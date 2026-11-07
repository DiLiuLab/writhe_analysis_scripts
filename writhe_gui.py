#!/usr/bin/env python3
r"""
Graphical front end shared by the writhe-analysis scripts.

Each analysis script opens this window when it is started without arguments
or with --gui:

    python analyze_projection_crossings.py            # empty form
    python analyze_projection_crossings.py --gui      # same
    python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz \
        --n-directions 1000 --gui                     # form pre-filled with these values

Running this module directly opens all four scripts as tabs of one window:

    python writhe_gui.py

The form is generated from each script's own argparse parser, so it offers
exactly the options documented by --help.  Run checks the values with that
parser, starts the script as a separate Python process, and streams its output
into the log.  The equivalent command line is shown above the log so a GUI run
can be repeated from a terminal.

Only the Python standard library is required.  tkinter is imported only when a
window is opened, so command-line use also works on Python builds without Tk.
"""

import argparse
import contextlib
import importlib
import io
import os
import queue
import shlex
import signal
import subprocess
import sys
import threading
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

GUI_FLAG = "--gui"

ParserFactory = Callable[[], argparse.ArgumentParser]
FormValue = Union[str, bool]
FileTypes = Tuple[Tuple[str, str], ...]

_XYZ_TYPES: FileTypes = (("XYZ coordinates", "*.xyz"), ("Text files", "*.txt"), ("All files", "*"))
_CSV_TYPES: FileTypes = (("CSV files", "*.csv"), ("All files", "*"))
_TXT_TYPES: FileTypes = (("Text files", "*.txt"), ("All files", "*"))

# Arguments that name files or folders (by argparse dest) and the dialog their
# Browse button opens: "open" an existing file, "save" a new file, or choose a
# "dir".  Add new path arguments here to give them a Browse button.
PATH_FIELDS: Dict[str, Tuple[str, FileTypes]] = {
    "xyz_file": ("open", _XYZ_TYPES),
    "xyz_a": ("open", _XYZ_TYPES),
    "xyz_b": ("open", _XYZ_TYPES),
    "scan_csv": ("open", _CSV_TYPES),
    "scan_a": ("open", _CSV_TYPES),
    "scan_b": ("open", _CSV_TYPES),
    "csv_out": ("save", _CSV_TYPES),
    "summary_out": ("save", _TXT_TYPES),
    "output": ("save", _TXT_TYPES),
    "outdir": ("dir", ()),
}

# The four scripts, in workflow order, for the tabbed window.
TOOLS = (
    ("1 Analyze", "analyze_projection_crossings"),
    ("2 Plot", "plot_writhe_projection_views"),
    ("3 Validate", "validate_projection_scan"),
    ("4 Compare", "compare_projection_scans"),
)

_COMMAND_COLOR = "#2f6fdf"
_OK_COLOR = "#2da44e"
_ERROR_COLOR = "#d1242f"
_NOTE_COLOR = "gray50"


# ---------------------------------------------------------------------------
# Command-line side: no tkinter needed.
# ---------------------------------------------------------------------------


def _asks_for_help(argv: Sequence[str]) -> bool:
    return any(arg == "-h" or (len(arg) >= 3 and "--help".startswith(arg)) for arg in argv)


def gui_requested(argv: Sequence[str]) -> bool:
    """
    True when a script should open the GUI: no arguments, or --gui among them.

    A request for help wins, so "--gui --help" prints the help text.  Scripts
    also open the GUI when argparse matches an abbreviation such as --gu.
    """
    if _asks_for_help(argv):
        return False
    return len(argv) == 0 or GUI_FLAG in argv


def add_gui_argument(parser: argparse.ArgumentParser) -> None:
    """Document --gui in a script's --help output."""
    parser.add_argument(
        GUI_FLAG,
        action="store_true",
        help="Open the graphical interface (also the default when no arguments are given); "
        "other arguments given with --gui pre-fill the form",
    )


@dataclass
class Field:
    """One command-line argument as it appears in the form."""

    dest: str
    option: Optional[str]  # long option string; None for a positional argument
    kind: str  # "text", "choice", "flag", or "multi" (one value per line)
    required: bool
    default: object
    help: str
    choices: Tuple[str, ...] = ()
    path_mode: Optional[str] = None  # "open", "save", "dir", or None
    file_types: FileTypes = ()
    const: object = None  # value stored by a flag

    @property
    def label(self) -> str:
        return self.option or self.dest


def _long_option(action: argparse.Action) -> Optional[str]:
    if not action.option_strings:
        return None
    for option in action.option_strings:
        if option.startswith("--"):
            return option
    return action.option_strings[0]


def _help_text(parser: argparse.ArgumentParser, action: argparse.Action) -> str:
    text = action.help or ""
    if "%" in text:
        params = dict(vars(action), prog=parser.prog)
        if params.get("choices") is not None:
            params["choices"] = ", ".join(str(choice) for choice in params["choices"])
        try:
            text = text % params
        except (KeyError, TypeError, ValueError):
            pass
    return " ".join(text.split())


def fields_from_parser(
    parser: argparse.ArgumentParser,
    path_fields: Optional[Dict[str, Tuple[str, FileTypes]]] = None,
) -> List[Field]:
    """Describe every user-settable argument of parser, in parser order."""
    path_fields = PATH_FIELDS if path_fields is None else path_fields
    fields: List[Field] = []
    for action in parser._actions:
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            continue
        if action.dest == "gui" or action.help == argparse.SUPPRESS:
            continue
        option = _long_option(action)
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            kind = "flag"
        elif type(action) is argparse._AppendAction and action.nargs is None:
            kind = "multi"
        elif type(action) is argparse._StoreAction and option is None and action.nargs in ("+", "*"):
            kind = "multi"
        elif type(action) is argparse._StoreAction and action.nargs is None:
            kind = "choice" if action.choices else "text"
        else:
            raise TypeError(
                "The GUI does not support argument {} ({}, nargs={!r})".format(
                    option or action.dest, type(action).__name__, action.nargs
                )
            )
        mode, file_types = path_fields.get(action.dest, (None, ()))
        fields.append(
            Field(
                dest=action.dest,
                option=option,
                kind=kind,
                required=bool(action.required),
                default=None if action.default is argparse.SUPPRESS else action.default,
                help=_help_text(parser, action),
                choices=tuple(str(choice) for choice in action.choices) if action.choices else (),
                path_mode=mode,
                file_types=tuple(file_types),
                const=action.const,
            )
        )
    return fields


def default_values(fields: Sequence[Field]) -> Dict[str, FormValue]:
    """Form values before the user changes anything."""
    values: Dict[str, FormValue] = {}
    for field in fields:
        if field.kind == "flag":
            values[field.dest] = False
        elif field.kind == "multi" or field.default is None:
            values[field.dest] = ""
        else:
            values[field.dest] = str(field.default)
    return values


def _option_tokens(option: str, value: str) -> List[str]:
    # "--direction=-0.27,0.86,0.42" keeps a value that starts with "-" from
    # being read as another option.
    if value.startswith("-"):
        return ["{}={}".format(option, value) if option.startswith("--") else option + value]
    return [option, value]


def build_argv(
    fields: Sequence[Field],
    values: Dict[str, FormValue],
    preview: bool = False,
) -> List[str]:
    """
    Translate form values into the argument list for the script.

    Blank optional values and values equal to their default are left out, so
    the command stays as short as the one a user would type.  A missing
    required value raises ValueError, or becomes a <name> placeholder when
    preview is True.
    """
    positionals: List[str] = []
    options: List[str] = []
    missing: List[str] = []
    for field in fields:
        value = values.get(field.dest)
        if field.kind == "flag":
            if value:
                options.append(str(field.option))
            continue
        if field.kind == "multi":
            items = [line.strip() for line in str(value or "").splitlines() if line.strip()]
        else:
            text = "" if value is None else str(value).strip()
            items = [text] if text else []
            if (
                items
                and field.option is not None
                and not field.required
                and field.default is not None
                and text == str(field.default)
            ):
                continue
        if not items:
            if field.required:
                missing.append(field.label)
                if preview and field.option is None:
                    positionals.append("<{}>".format(field.label))
                elif preview:
                    options.extend([field.option, "<value>"])
            continue
        if field.option is None:
            positionals.extend(items)
        else:
            for item in items:
                options.extend(_option_tokens(field.option, item))
    if missing and not preview:
        raise ValueError("Please fill in: {}".format(", ".join(missing)))
    if any(item.startswith("-") for item in positionals):
        # "--" ends the options, so a file name such as "-curve.xyz" stays a positional.
        return options + ["--"] + positionals
    return positionals + options


def _last_error_line(captured: str) -> str:
    lines = [line.strip() for line in captured.splitlines() if line.strip()]
    if not lines:
        return "invalid arguments"
    message = lines[-1]
    marker = ": error: "
    return message.split(marker, 1)[1] if marker in message else message


def check_argv(parser: argparse.ArgumentParser, argv: Sequence[str]) -> Optional[str]:
    """Return argparse's error message for argv, or None when it parses."""
    parser.color = False  # Python 3.15+ would otherwise add ANSI color codes
    captured = io.StringIO()
    try:
        with contextlib.redirect_stderr(captured), contextlib.redirect_stdout(captured):
            parser.parse_args(list(argv))
    except SystemExit as exc:
        if not exc.code:
            return None
        return _last_error_line(captured.getvalue())
    return None


def initial_values(
    parser_factory: ParserFactory,
    argv: Sequence[str],
) -> Tuple[Dict[str, FormValue], Optional[str]]:
    """
    Form values: the defaults, overridden by any arguments given with --gui.

    Returns the values and, when the arguments could not be read, a warning.
    Values are kept as typed (for example, a --direction vector is not
    normalized) and required arguments may be missing.
    """
    fields = fields_from_parser(parser_factory())
    values = default_values(fields)
    given = [arg for arg in argv if arg != GUI_FLAG]
    if not given:
        return values, None

    lenient = parser_factory()
    for action in lenient._actions:
        action.required = False
        if not action.option_strings:
            if action.nargs is None:
                action.nargs = "?"
            elif action.nargs == "+":
                action.nargs = "*"
        if type(action) in (argparse._StoreAction, argparse._AppendAction):
            action.type = None
    lenient.color = False
    namespace = None
    error = ""
    # parse_intermixed_args accepts positionals between options; on Python 3.9
    # and 3.12.0 it rejects a "--" followed by a value starting with "-", which
    # plain parse_args handles.
    for parse in (lenient.parse_intermixed_args, lenient.parse_args):
        captured = io.StringIO()
        try:
            with contextlib.redirect_stderr(captured), contextlib.redirect_stdout(captured):
                namespace = parse(given)
            break
        except SystemExit as exc:
            if not exc.code:  # -h/--help: nothing to pre-fill
                return values, None
            error = _last_error_line(captured.getvalue())
    if namespace is None:
        return values, "Ignored the command-line arguments ({}).".format(error)

    for field in fields:
        raw = getattr(namespace, field.dest, None)
        if field.kind == "flag":
            values[field.dest] = raw == field.const
        elif field.kind == "multi":
            if raw:
                values[field.dest] = "\n".join(str(item) for item in raw)
        elif raw is not None:
            values[field.dest] = str(raw)
    return values, None


def relative_to(path: str, base: str) -> str:
    """path relative to base when it lies inside base; otherwise absolute."""
    absolute = os.path.abspath(path)
    base = os.path.abspath(base)
    if os.path.dirname(base) == base:
        return absolute
    try:
        relative = os.path.relpath(absolute, base)
    except ValueError:  # another drive on Windows
        return absolute
    if relative == os.pardir or relative.startswith(os.pardir + os.sep):
        return absolute
    return relative


_CMD_SPECIAL = frozenset(" \t&|<>^()%!\"")


def quote_windows(arg: str) -> str:
    """
    Quote one argument for cmd.exe and the Microsoft C runtime.

    Unlike subprocess.list2cmdline, this also quotes characters that cmd.exe
    treats specially, such as '&' (legal in Windows file names).
    """
    if arg and not any(ch in _CMD_SPECIAL for ch in arg):
        return arg
    out, backslashes = ['"'], 0
    for ch in arg:
        if ch == "\\":
            backslashes += 1
            continue
        if ch == '"':
            out.append("\\" * (2 * backslashes + 1) + '"')
        else:
            out.append("\\" * backslashes + ch)
        backslashes = 0
    out.append("\\" * (2 * backslashes) + '"')
    return "".join(out)


def format_command(script_path: str, argv: Sequence[str], cwd: Optional[str] = None) -> str:
    """The command line a user could type in cwd to repeat a run."""
    words = ["python", relative_to(script_path, cwd or os.getcwd())] + list(argv)
    if os.name == "nt":
        return " ".join(quote_windows(word) for word in words)
    return shlex.join(words)


def _pump_output(process: "subprocess.Popen[str]", sink: "queue.Queue[Tuple[str, object]]") -> None:
    """Background thread: forward the child's output lines, then its exit code."""
    try:
        assert process.stdout is not None
        for line in process.stdout:
            sink.put(("output", line))
    except (OSError, ValueError):
        pass
    finally:
        try:
            if process.stdout is not None:
                process.stdout.close()
        except OSError:
            pass
        sink.put(("exit", process.wait()))


# ---------------------------------------------------------------------------
# Window side: tkinter is imported here, only when a window is opened.
# ---------------------------------------------------------------------------


def _tk():
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    return tk, ttk, filedialog, messagebox


def _create_root(title: str):
    """Return (Tk root, None), or (None, reason) when no window can be opened."""
    try:
        import tkinter as tk
    except ImportError as exc:
        return None, "tkinter is not available in this Python installation ({})".format(exc)
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        return None, "no graphical display is available ({})".format(exc)
    root.title(title)
    return root, None


class ScriptForm:
    """Parameter form, command preview, Run/Stop buttons, and log for one script."""

    POLL_MS = 100

    def __init__(
        self,
        master,
        parser_factory: ParserFactory,
        script_path: str,
        argv: Sequence[str] = (),
        path_fields: Optional[Dict[str, Tuple[str, FileTypes]]] = None,
        workdir: Optional[str] = None,
    ) -> None:
        tk, ttk, filedialog, messagebox = _tk()
        self._tk = tk
        self._ttk = ttk
        self._filedialog = filedialog
        self._messagebox = messagebox
        self.parser_factory = parser_factory
        self.script_path = os.path.abspath(script_path)
        self.workdir = os.path.abspath(workdir or os.getcwd())
        parser = parser_factory()
        self.fields = fields_from_parser(parser, path_fields)
        values, warning = initial_values(parser_factory, list(argv))

        self.process: Optional[subprocess.Popen] = None
        self.returncode: Optional[int] = None
        self._stopping = False
        self._queue: "queue.Queue[Tuple[str, object]]" = queue.Queue()
        self._poll_job = None
        self._vars: Dict[str, object] = {}
        self._texts: Dict[str, object] = {}
        self._field_labels: List[object] = []
        self._help_labels: List[object] = []

        style = ttk.Style(master)
        style.configure("Help.TLabel", foreground=_NOTE_COLOR, font="TkSmallCaptionFont")

        self.frame = ttk.Frame(master, padding=10)
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(2, weight=1)
        ttk.Label(self.frame, text=os.path.basename(self.script_path), font="TkHeadingFont").grid(
            row=0, column=0, sticky="w"
        )
        header = (parser.description or "").strip()
        header += "\nRelative paths are resolved from the working folder: {}".format(self.workdir)
        header_label = ttk.Label(self.frame, text=header, wraplength=820, justify="left")
        header_label.grid(row=1, column=0, sticky="w", pady=(2, 0))
        # Rewrap the header when the window is resized instead of cutting it off.
        self.frame.bind(
            "<Configure>",
            lambda event: header_label.configure(wraplength=max(200, event.width - 24)),
            add="+",
        )
        paned = ttk.Panedwindow(self.frame, orient="vertical")
        paned.grid(row=2, column=0, sticky="nsew", pady=(8, 0))
        paned.add(self._build_fields(paned, values), weight=3)
        paned.add(self._build_runner(paned), weight=2)
        paned.bind("<Configure>", lambda event: self._keep_runner_visible(paned, event.height), add="+")

        self._log(
            "Python: {}\nWorking folder: {}\n".format(sys.executable, self.workdir), "note"
        )
        if warning:
            self._log(warning + "\n", "error")
        self.update_preview()

    # -- layout -------------------------------------------------------------

    def _build_fields(self, parent, values: Dict[str, FormValue]):
        tk, ttk = self._tk, self._ttk
        outer = ttk.Frame(parent)
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(0, weight=1)
        canvas = tk.Canvas(outer, highlightthickness=0, borderwidth=0)
        background = ttk.Style(parent).lookup("TFrame", "background")
        if background:
            try:
                canvas.configure(background=background)
            except tk.TclError:
                pass
        scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        inner = ttk.Frame(canvas, padding=(2, 0, 10, 6))
        window = canvas.create_window(0, 0, window=inner, anchor="nw")
        inner.bind("<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all")))

        def fit_width(event) -> None:
            canvas.itemconfigure(window, width=event.width)
            label_width = max((label.winfo_reqwidth() for label in self._field_labels), default=0)
            wrap = max(150, event.width - label_width - 40)
            for help_label in self._help_labels:
                help_label.configure(wraplength=wrap)

        canvas.bind("<Configure>", fit_width)
        inner.columnconfigure(1, weight=1)

        row = 0
        for field in self.fields:
            row = self._add_field(inner, row, field, values[field.dest])
        ttk.Label(inner, text="* required", style="Help.TLabel").grid(
            row=row, column=1, sticky="w", pady=(8, 0)
        )
        self._enable_wheel(canvas)
        inner.update_idletasks()
        canvas.configure(
            width=inner.winfo_reqwidth(),
            height=min(inner.winfo_reqheight(), int(0.48 * canvas.winfo_screenheight())),
        )
        return outer

    def _add_field(self, parent, row: int, field: Field, value: FormValue) -> int:
        tk, ttk = self._tk, self._ttk
        label = ttk.Label(parent, text=field.label + (" *" if field.required else ""), font="TkFixedFont")
        label.grid(row=row, column=0, sticky="nw", padx=(0, 12), pady=(8, 0))
        self._field_labels.append(label)
        var = None
        if field.kind == "flag":
            var = tk.BooleanVar(master=parent, value=bool(value))
            ttk.Checkbutton(parent, variable=var).grid(row=row, column=1, sticky="w", pady=(8, 0))
        elif field.kind == "choice":
            var = tk.StringVar(master=parent, value=str(value))
            ttk.Combobox(
                parent, textvariable=var, values=list(field.choices), state="readonly", width=12
            ).grid(row=row, column=1, sticky="w", pady=(8, 0))
        elif field.kind == "multi":
            text = tk.Text(parent, height=3, width=52, wrap="none", undo=True, font="TkFixedFont")
            text.insert("1.0", str(value))
            text.edit_modified(False)
            text.bind("<<Modified>>", self._on_text_modified)
            text.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(8, 0))
            self._texts[field.dest] = text
        else:
            var = tk.StringVar(master=parent, value=str(value))
            entry = ttk.Entry(parent, textvariable=var, width=52)
            if field.path_mode:
                entry.grid(row=row, column=1, sticky="ew", pady=(8, 0))
                ttk.Button(parent, text="Browse…", command=lambda f=field: self._browse(f)).grid(
                    row=row, column=2, sticky="w", padx=(6, 0), pady=(8, 0)
                )
            else:
                entry.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(8, 0))
        if var is not None:
            var.trace_add("write", lambda *_args: self.update_preview())
            self._vars[field.dest] = var

        help_text = field.help
        if field.kind == "multi":
            help_text = "One value per line. " + help_text
        if help_text:
            help_label = ttk.Label(parent, text=help_text, style="Help.TLabel", wraplength=600, justify="left")
            help_label.grid(row=row + 1, column=1, columnspan=2, sticky="w")
            self._help_labels.append(help_label)
        return row + 2

    def _build_runner(self, parent):
        tk, ttk = self._tk, self._ttk
        box = ttk.Frame(parent)
        box.columnconfigure(0, weight=1)
        box.rowconfigure(3, weight=1)
        ttk.Label(box, text="Equivalent command").grid(row=0, column=0, sticky="w", pady=(8, 2))
        self.command_text = tk.Text(box, height=3, wrap="word", font="TkFixedFont", takefocus=0)
        self.command_text.grid(row=1, column=0, sticky="ew")

        buttons = ttk.Frame(box)
        buttons.grid(row=2, column=0, sticky="ew", pady=6)
        self.run_button = ttk.Button(buttons, text="Run", command=self.run)
        self.stop_button = ttk.Button(buttons, text="Stop", command=self.stop, state="disabled")
        self.run_button.pack(side="left")
        self.stop_button.pack(side="left", padx=(6, 0))
        ttk.Button(buttons, text="Copy command", command=self.copy_command).pack(
            side="left", padx=(18, 0)
        )
        ttk.Button(buttons, text="Clear log", command=self.clear_log).pack(side="left", padx=(6, 0))
        self.status = ttk.Label(buttons, text="Ready")
        self.status.pack(side="left", padx=(18, 0))

        log_frame = ttk.Frame(box)
        log_frame.grid(row=3, column=0, sticky="nsew")
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, height=12, wrap="word", font="TkFixedFont", state="disabled")
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.log_text.tag_configure("command", foreground=_COMMAND_COLOR)
        self.log_text.tag_configure("ok", foreground=_OK_COLOR)
        self.log_text.tag_configure("error", foreground=_ERROR_COLOR)
        self.log_text.tag_configure("note", foreground=_NOTE_COLOR)
        return box

    def _keep_runner_visible(self, paned, height: int, smallest: int = 210) -> None:
        """When the window gets shorter, shrink the form rather than the Run buttons and log."""
        limit = height - smallest
        try:
            if limit > 0 and paned.sashpos(0) > limit:
                paned.sashpos(0, limit)
        except self._tk.TclError:
            pass

    def _enable_wheel(self, canvas) -> None:
        """Scroll the form with the mouse wheel while the pointer is over it."""
        tk = self._tk

        def scroll(event):
            widget = event.widget
            node = widget
            while isinstance(node, tk.Misc) and node is not canvas:
                node = node.master
            if node is not canvas or canvas.yview() == (0.0, 1.0):
                return
            if isinstance(widget, tk.Text) and widget.yview() != (0.0, 1.0):
                return  # a multi-line box with more lines than it shows scrolls itself
            if getattr(event, "num", None) == 4:
                step = -1
            elif getattr(event, "num", None) == 5:
                step = 1
            elif sys.platform == "darwin":
                step = -event.delta
            else:
                step = -int(event.delta / 120) or (-1 if event.delta > 0 else 1)
            canvas.yview_scroll(step, "units")

        for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            canvas.bind_all(sequence, scroll, add="+")
        # Tk's combobox steps through its values under the wheel, so scrolling
        # the form across --engine or --format would silently change them.
        for sequence in canvas.bind_class("TCombobox"):
            if "MouseWheel" in sequence or sequence in ("<Button-4>", "<Button-5>"):
                canvas.bind_class("TCombobox", sequence, "")

    # -- values -------------------------------------------------------------

    def _collect_values(self) -> Dict[str, FormValue]:
        values: Dict[str, FormValue] = {}
        for field in self.fields:
            if field.dest in self._texts:
                values[field.dest] = self._texts[field.dest].get("1.0", "end-1c")
                continue
            value = self._vars[field.dest].get()
            if field.path_mode and isinstance(value, str) and value.strip().startswith("~"):
                value = os.path.expanduser(value.strip())  # no shell expands "~" for us
            values[field.dest] = value
        return values

    def _resolve(self, path: str) -> str:
        return os.path.join(self.workdir, os.path.expanduser(path))

    def _missing_inputs(self, values: Dict[str, FormValue]) -> Optional[str]:
        for field in self.fields:
            if field.path_mode != "open":
                continue
            text = str(values.get(field.dest) or "").strip()
            if text and not os.path.isfile(self._resolve(text)):
                return "File not found for {}:\n{}".format(field.label, text)
        return None

    def _browse(self, field: Field) -> None:
        current = str(self._vars[field.dest].get()).strip()
        start = self._resolve(current) if current else self.workdir
        initialdir = start if os.path.isdir(start) else os.path.dirname(start)
        if not os.path.isdir(initialdir):
            initialdir = self.workdir
        options = {
            "parent": self.frame.winfo_toplevel(),
            "initialdir": initialdir,
            "title": "Choose {}".format(field.label),
        }
        if field.path_mode == "dir":
            chosen = self._filedialog.askdirectory(mustexist=False, **options)
        elif field.path_mode == "save":
            pattern = field.file_types[0][1] if field.file_types else ""
            extension = pattern[1:] if pattern.startswith("*.") else ""
            chosen = self._filedialog.asksaveasfilename(
                filetypes=list(field.file_types),
                defaultextension=extension,
                initialfile=os.path.basename(current),
                **options,
            )
        else:
            chosen = self._filedialog.askopenfilename(filetypes=list(field.file_types), **options)
        if chosen:
            self._vars[field.dest].set(relative_to(chosen, self.workdir))

    def _on_text_modified(self, event) -> None:
        if event.widget.edit_modified():
            event.widget.edit_modified(False)
            self.update_preview()

    def update_preview(self) -> None:
        argv = build_argv(self.fields, self._collect_values(), preview=True)
        self.command_text.configure(state="normal")
        self.command_text.delete("1.0", "end")
        self.command_text.insert("1.0", format_command(self.script_path, argv, self.workdir))
        self.command_text.configure(state="disabled")

    # -- running ------------------------------------------------------------

    @property
    def running(self) -> bool:
        return self.process is not None

    def run(self) -> None:
        """Validate the form and start the script."""
        if self.running:
            return
        values = self._collect_values()
        try:
            argv = build_argv(self.fields, values)
        except ValueError as exc:
            self._messagebox.showerror("Missing values", str(exc), parent=self.frame.winfo_toplevel())
            return
        problem = check_argv(self.parser_factory(), argv) or self._missing_inputs(values)
        if problem:
            self._messagebox.showerror("Cannot run", problem, parent=self.frame.winfo_toplevel())
            return
        self.start(argv)

    def start(self, argv: Sequence[str]) -> None:
        """Run the script with argv in a child process (no validation)."""
        if self.running:
            return
        command = [sys.executable, "-u", self.script_path] + list(argv)
        env = dict(os.environ)
        # The scripts only save figures; Agg avoids opening any plotting window.
        env.update({
            "MPLBACKEND": "Agg",
            "PYTHONUNBUFFERED": "1",
            "PYTHONIOENCODING": "utf-8",
            "PYTHON_COLORS": "0",  # plain error text in the log on Python 3.13+
        })
        self._log("\n$ {}\n".format(format_command(self.script_path, argv, self.workdir)), "command")
        try:
            self.process = subprocess.Popen(
                command,
                cwd=self.workdir,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=env,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
        except OSError as exc:
            self.process = None
            self._log("Could not start the script: {}\n".format(exc), "error")
            return
        self.returncode = None
        self._stopping = False
        self._queue = queue.Queue()
        threading.Thread(target=_pump_output, args=(self.process, self._queue), daemon=True).start()
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.status.configure(text="Running…")
        self._poll_job = self.frame.after(self.POLL_MS, self._poll)

    def _poll(self) -> None:
        self._poll_job = None
        chunks: List[str] = []
        exit_code = None
        try:
            for _ in range(2000):
                kind, payload = self._queue.get_nowait()
                if kind == "exit":
                    exit_code = payload
                    break
                chunks.append(str(payload))
        except queue.Empty:
            pass
        if chunks:
            self._log("".join(chunks))
        if exit_code is None:
            self._poll_job = self.frame.after(self.POLL_MS, self._poll)
        else:
            self._finish(int(exit_code))

    def _finish(self, code: int) -> None:
        self.process = None
        self.returncode = code
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        if code == 0:
            message, tag = "Finished successfully.", "ok"
        elif self._stopping:
            message, tag = "Stopped.", "error"
        elif code < 0:  # killed by a signal the GUI did not send (crash, out of memory, ...)
            try:
                name = signal.Signals(-code).name
            except ValueError:
                name = "signal {}".format(-code)
            message, tag = "Terminated by {} (exit code {}); see the messages above.".format(name, code), "error"
        else:
            message, tag = "Failed with exit code {}; see the messages above.".format(code), "error"
        self._log(message + "\n", tag)
        self.status.configure(text=message)

    def stop(self) -> None:
        """Ask the running script to stop; force it after three seconds."""
        process = self.process
        if process is None or process.poll() is not None:
            return
        self._stopping = True
        self._log("Stopping…\n", "error")
        try:
            process.terminate()
        except OSError:
            pass
        self.frame.after(3000, lambda: self._kill(process))

    @staticmethod
    def _kill(process: subprocess.Popen) -> None:
        if process.poll() is None:
            try:
                process.kill()
            except OSError:
                pass

    def shutdown(self) -> None:
        """Kill any running script and stop polling (used when the window closes)."""
        process, self.process = self.process, None
        if self._poll_job is not None:
            try:
                self.frame.after_cancel(self._poll_job)
            except Exception:
                pass
            self._poll_job = None
        if process is not None and process.poll() is None:
            self._kill(process)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass

    # -- small actions ------------------------------------------------------

    def copy_command(self) -> None:
        top = self.frame.winfo_toplevel()
        top.clipboard_clear()
        top.clipboard_append(self.command_text.get("1.0", "end-1c"))
        self.status.configure(text="Command copied to the clipboard.")

    def clear_log(self) -> None:
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

    def _log(self, text: str, tag: Optional[str] = None) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", text, (tag,) if tag else ())
        self.log_text.see("end")
        self.log_text.configure(state="disabled")


def _fit_to_screen(root) -> None:
    root.update_idletasks()
    width = min(max(root.winfo_reqwidth(), 760), int(0.95 * root.winfo_screenwidth()))
    height = min(root.winfo_reqheight(), int(0.88 * root.winfo_screenheight()))
    root.geometry("{}x{}".format(width, height))
    root.minsize(640, 460)


def _run_window(root, forms: Sequence[ScriptForm]) -> int:
    """Show root until it is closed; stop any running script on the way out."""
    _tk_module, _ttk, _filedialog, messagebox = _tk()

    def close() -> None:
        if any(form.running for form in forms) and not messagebox.askyesno(
            "Quit", "A run is still in progress. Stop it and quit?", parent=root
        ):
            return
        for form in forms:
            form.shutdown()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", close)
    if sys.platform == "darwin":
        root.createcommand("tk::mac::Quit", close)

    # Tk installs its own SIGINT/SIGTERM/SIGHUP handlers on macOS, which end the
    # process without stopping a running script.  Stop the scripts first, then
    # leave the main loop.  (Raising KeyboardInterrupt would not work: Tkinter
    # swallows exceptions raised inside its callbacks.)
    received: List[int] = []

    def on_signal(signum, _frame) -> None:
        received.append(signum)
        for form in forms:
            form.shutdown()
        root.quit()

    previous = {}
    for name in ("SIGINT", "SIGTERM", "SIGHUP"):
        number = getattr(signal, name, None)
        if number is None:
            continue
        try:
            current = signal.getsignal(number)
            signal.signal(number, signal.SIG_IGN if current == signal.SIG_IGN else on_signal)
            previous[number] = current
        except (ValueError, OSError):  # not the main thread
            pass

    _fit_to_screen(root)
    try:
        root.mainloop()
    except KeyboardInterrupt:
        pass
    finally:
        for form in forms:
            form.shutdown()
        for number, handler in previous.items():
            try:
                signal.signal(number, handler if handler is not None else signal.SIG_DFL)
            except (ValueError, OSError, TypeError):
                pass
        if received:
            try:
                root.destroy()
            except Exception:
                pass
    return 128 + received[0] if received else 0


def launch_gui(parser_factory: ParserFactory, script_path: str, argv: Sequence[str] = ()) -> int:
    """
    Open the parameter window for one script and return an exit status.

    argv may hold the script's own arguments (with or without --gui); they
    pre-fill the form.  Without a usable display, print the reason and the
    script's --help text and return 2.
    """
    root, problem = _create_root("Writhe analysis: {}".format(os.path.basename(script_path)))
    if root is None:
        sys.stderr.write(
            "Cannot open the graphical interface: {}.\n"
            "Run the script with command-line arguments instead:\n\n{}".format(
                problem, parser_factory().format_help()
            )
        )
        return 2
    form = ScriptForm(root, parser_factory, script_path, argv)
    form.frame.pack(fill="both", expand=True)
    return _run_window(root, [form])


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Open all four scripts as tabs of one window."""
    argv = [arg for arg in (sys.argv[1:] if argv is None else argv) if arg != GUI_FLAG]
    if argv:
        print((__doc__ or "").strip())
        return 0 if argv[0] in ("-h", "--help") else 2
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    os.environ.setdefault("MPLBACKEND", "Agg")
    root, problem = _create_root("Writhe analysis scripts")
    if root is None:
        sys.stderr.write(
            "Cannot open the graphical interface: {}.\n"
            "Run the individual scripts with --help to see their command-line options.\n".format(problem)
        )
        return 2
    _tk_module, ttk, _filedialog, _messagebox = _tk()
    notebook = ttk.Notebook(root)
    notebook.pack(fill="both", expand=True)
    forms: List[ScriptForm] = []
    for title, module_name in TOOLS:
        try:
            module = importlib.import_module(module_name)
            form = ScriptForm(notebook, module.build_parser, module.__file__)
        except Exception as exc:  # e.g. a missing third-party package
            frame = ttk.Frame(notebook, padding=20)
            ttk.Label(
                frame, text="Cannot load {}.py:\n{}".format(module_name, exc), wraplength=700
            ).pack(anchor="w")
            notebook.add(frame, text=title)
            continue
        notebook.add(form.frame, text=title)
        forms.append(form)
    return _run_window(root, forms)


if __name__ == "__main__":
    sys.exit(main())
