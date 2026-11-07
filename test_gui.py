#!/usr/bin/env python3
"""
Tests for the graphical interface in writhe_gui.py.

The argument-translation and launch tests need no display.  The window tests
run only when Tk can open a window; the window stays hidden.
"""

import contextlib
import io
import os
import shlex
import tempfile
import time
import unittest
from unittest import mock

import numpy as np

import analyze_projection_crossings as analyze
import compare_projection_scans as compare
import plot_writhe_projection_views as plot
import validate_projection_scan as validate
import writhe_gui as gui


HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "examples")
XYZ = os.path.join(EXAMPLES, "sc_W3_ALocal.xyz")
SCRIPTS = (analyze, plot, validate, compare)


def form_values(module, **changes):
    """Default form values for a script with some fields changed."""
    fields = gui.fields_from_parser(module.build_parser())
    values = gui.default_values(fields)
    values.update(changes)
    return fields, values


def roundtrip(module, **changes):
    """Build argv from form values and parse it with the script's own parser."""
    fields, values = form_values(module, **changes)
    argv = gui.build_argv(fields, values)
    return argv, module.build_parser().parse_args(argv)


class FieldTests(unittest.TestCase):
    def test_every_command_line_option_has_a_field(self) -> None:
        for module in SCRIPTS:
            parser = module.build_parser()
            expected = {action.dest for action in parser._actions} - {"help", "gui"}
            found = {field.dest for field in gui.fields_from_parser(parser)}
            self.assertEqual(found, expected, module.__name__)

    def test_field_kinds_and_dialogs(self) -> None:
        fields = {field.dest: field for field in gui.fields_from_parser(analyze.build_parser())}
        self.assertEqual((fields["xyz_file"].kind, fields["xyz_file"].path_mode), ("text", "open"))
        self.assertTrue(fields["xyz_file"].required)
        self.assertEqual(fields["include_axis_views"].kind, "flag")
        self.assertEqual(fields["engine"].kind, "choice")
        self.assertEqual(fields["engine"].choices, ("auto", "shapely", "numpy"))
        self.assertEqual(fields["direction"].kind, "multi")
        self.assertEqual(fields["csv_out"].path_mode, "save")
        plot_fields = {field.dest: field for field in gui.fields_from_parser(plot.build_parser())}
        self.assertEqual(plot_fields["outdir"].path_mode, "dir")

    def test_unsupported_action_is_reported(self) -> None:
        parser = analyze.build_parser()
        parser.add_argument("--verbose", "-v", action="count")
        with self.assertRaises(TypeError):
            gui.fields_from_parser(parser)


class ArgvTests(unittest.TestCase):
    def test_defaults_are_left_out(self) -> None:
        argv, args = roundtrip(analyze, xyz_file="examples/sc_W3_ALocal.xyz")
        self.assertEqual(argv, ["examples/sc_W3_ALocal.xyz"])
        self.assertEqual(args.n_directions, 10000)
        self.assertEqual(args.engine, "auto")

    def test_changed_values_flags_choices_and_zero(self) -> None:
        argv, args = roundtrip(
            analyze,
            xyz_file="a.xyz",
            n_directions=" 1000 ",
            analysis_points="0",
            engine="numpy",
            include_axis_views=True,
        )
        self.assertEqual(
            argv,
            ["a.xyz", "--n-directions", "1000", "--analysis-points", "0", "--engine", "numpy", "--include-axis-views"],
        )
        self.assertEqual((args.n_directions, args.analysis_points, args.engine), (1000, 0, "numpy"))
        self.assertTrue(args.include_axis_views)

    def test_negative_vectors_and_repeated_options(self) -> None:
        argv, args = roundtrip(
            analyze,
            xyz_file="a.xyz",
            direction="0.2494,0.8564,0.4521\n\n  -0.27,0.8643,0.4243  \n",
        )
        self.assertIn("--direction=-0.27,0.8643,0.4243", argv)
        self.assertEqual(len(args.direction), 2)
        self.assertLess(args.direction[1][0], 0.0)
        np.testing.assert_allclose(np.linalg.norm(args.direction[1]), 1.0)

    def test_named_views(self) -> None:
        _argv, args = roundtrip(plot, xyz_file="a.xyz", view="my_view:-0.2,0.3,0.9\nother:0,0,1")
        self.assertEqual([name for name, _normal in args.view], ["my_view", "other"])
        self.assertLess(args.view[0][1][0], 0.0)

    def test_compare_positionals_keep_their_order(self) -> None:
        argv, args = roundtrip(
            compare,
            xyz_a="a.xyz",
            scan_a="a.csv",
            xyz_b="b.xyz",
            scan_b="b.csv",
            label_a="ALocal",
            label_b="R11",
            output_format="both",
        )
        self.assertEqual(argv[:4], ["a.xyz", "a.csv", "b.xyz", "b.csv"])
        self.assertEqual((args.xyz_a, args.scan_a, args.xyz_b, args.scan_b), ("a.xyz", "a.csv", "b.xyz", "b.csv"))
        self.assertEqual((args.label_a, args.label_b, args.output_format), ("ALocal", "R11", "both"))

    def test_paths_with_spaces_survive_the_displayed_command(self) -> None:
        path = os.path.join("My Data", "curve one.xyz")
        argv, args = roundtrip(validate, xyz_file=path, scan_csv="scan 1.csv", output="report (v2).txt")
        self.assertEqual((args.xyz_file, args.scan_csv, args.output), (path, "scan 1.csv", "report (v2).txt"))
        command = gui.format_command(os.path.join(HERE, "validate_projection_scan.py"), argv, HERE)
        if os.name != "nt":
            self.assertEqual(shlex.split(command), ["python", "validate_projection_scan.py"] + argv)

    def test_file_name_starting_with_dash(self) -> None:
        argv, args = roundtrip(analyze, xyz_file="-odd.xyz", n_directions="10")
        self.assertIn("--", argv)
        self.assertEqual((args.xyz_file, args.n_directions), ("-odd.xyz", 10))

    def test_missing_required_value(self) -> None:
        fields, values = form_values(compare, xyz_a="a.xyz")
        with self.assertRaises(ValueError) as caught:
            gui.build_argv(fields, values)
        self.assertIn("scan_a", str(caught.exception))
        preview = gui.build_argv(fields, values, preview=True)
        self.assertEqual(preview[:4], ["a.xyz", "<scan_a>", "<xyz_b>", "<scan_b>"])

    def test_check_argv_reports_parser_errors(self) -> None:
        fields, values = form_values(analyze, xyz_file="a.xyz", n_directions="ten")
        message = gui.check_argv(analyze.build_parser(), gui.build_argv(fields, values))
        self.assertIsNotNone(message)
        self.assertIn("invalid int value", message)
        fields, values = form_values(analyze, xyz_file="a.xyz", direction="1,2")
        self.assertIsNotNone(gui.check_argv(analyze.build_parser(), gui.build_argv(fields, values)))
        self.assertIsNone(gui.check_argv(analyze.build_parser(), ["a.xyz"]))

    def test_relative_to(self) -> None:
        self.assertEqual(gui.relative_to(XYZ, HERE), os.path.join("examples", "sc_W3_ALocal.xyz"))
        self.assertEqual(gui.relative_to(HERE, EXAMPLES), HERE)


class LaunchTests(unittest.TestCase):
    def test_gui_requested(self) -> None:
        self.assertTrue(gui.gui_requested([]))
        self.assertTrue(gui.gui_requested(["--gui"]))
        self.assertTrue(gui.gui_requested(["a.xyz", "--gui"]))
        self.assertFalse(gui.gui_requested(["a.xyz"]))
        self.assertFalse(gui.gui_requested(["--help"]))

    def test_gui_option_is_documented(self) -> None:
        for module in SCRIPTS:
            self.assertIn("--gui", module.build_parser().format_help(), module.__name__)

    def test_scripts_open_the_gui_without_arguments_or_with_gui(self) -> None:
        for module in SCRIPTS:
            for argv in ([], ["--gui"], ["examples/sc_W3_ALocal.xyz", "--gui"]):
                with mock.patch.object(module, "launch_gui", return_value=0) as launch:
                    with self.assertRaises(SystemExit) as caught:
                        module.main(argv)
                self.assertEqual(caught.exception.code, 0)
                launch.assert_called_once()
                factory, script, passed = launch.call_args[0]
                self.assertIs(factory, module.build_parser)
                self.assertEqual(os.path.abspath(script), os.path.abspath(module.__file__))
                self.assertEqual(list(passed), argv)

    def test_command_line_use_does_not_open_the_gui(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "scan.csv")
            summary_path = os.path.join(tmp, "summary.txt")
            with mock.patch.object(analyze, "launch_gui") as launch, contextlib.redirect_stdout(io.StringIO()):
                analyze.main([
                    XYZ, "--n-directions", "20", "--analysis-points", "200",
                    "--progress-every", "0", "--csv-out", csv_path, "--summary-out", summary_path,
                ])
            launch.assert_not_called()
            self.assertTrue(os.path.isfile(csv_path))
            self.assertTrue(os.path.isfile(summary_path))

    def test_prefill_from_command_line(self) -> None:
        values, warning = gui.initial_values(
            analyze.build_parser,
            [
                "examples/sc_W3_ALocal.xyz", "--n-directions", "1000",
                "--direction=-0.27,0.86,0.42", "--direction", "0.1,0.2,0.3",
                "--include-axis-views", "--gui",
            ],
        )
        self.assertIsNone(warning)
        self.assertEqual(values["xyz_file"], "examples/sc_W3_ALocal.xyz")
        self.assertEqual(values["n_directions"], "1000")
        self.assertEqual(values["analysis_points"], "1000")
        self.assertEqual(values["direction"], "-0.27,0.86,0.42\n0.1,0.2,0.3")
        self.assertIs(values["include_axis_views"], True)
        self.assertEqual(values["engine"], "auto")

    def test_prefill_with_positionals_between_options(self) -> None:
        values, warning = gui.initial_values(
            compare.build_parser, ["--label-a", "X", "a.xyz", "--gui", "a.csv"]
        )
        self.assertIsNone(warning)
        self.assertEqual(
            (values["xyz_a"], values["scan_a"], values["xyz_b"], values["label_a"]),
            ("a.xyz", "a.csv", "", "X"),
        )

    def test_prefill_ignores_unreadable_arguments(self) -> None:
        values, warning = gui.initial_values(analyze.build_parser, ["--bogus", "--gui"])
        self.assertIsNotNone(warning)
        self.assertEqual(values, gui.default_values(gui.fields_from_parser(analyze.build_parser())))

    def test_without_a_display_the_help_is_printed(self) -> None:
        stderr = io.StringIO()
        with mock.patch.object(gui, "_create_root", return_value=(None, "no display")), contextlib.redirect_stderr(stderr):
            code = gui.launch_gui(analyze.build_parser, analyze.__file__, [])
        self.assertEqual(code, 2)
        self.assertIn("no display", stderr.getvalue())
        self.assertIn("--n-directions", stderr.getvalue())


class WindowTests(unittest.TestCase):
    """Build the real forms in a hidden window and run one short analysis."""

    @classmethod
    def setUpClass(cls) -> None:
        try:
            import tkinter
        except ImportError:
            raise unittest.SkipTest("tkinter is not available")
        try:
            cls.root = tkinter.Tk()
        except tkinter.TclError as exc:
            raise unittest.SkipTest("no display: {}".format(exc))
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.root.destroy()

    def make_form(self, module, argv=()):
        form = gui.ScriptForm(self.root, module.build_parser, module.__file__, argv, workdir=HERE)
        self.addCleanup(form.frame.destroy)
        self.addCleanup(form.shutdown)
        return form

    def preview(self, form) -> str:
        return form.command_text.get("1.0", "end-1c")

    def test_every_script_builds_a_form(self) -> None:
        for module in SCRIPTS:
            form = self.make_form(module)
            self.assertIn(os.path.basename(module.__file__), self.preview(form))

    def test_preview_follows_the_fields(self) -> None:
        form = self.make_form(analyze, ["examples/sc_W3_ALocal.xyz", "--gui"])
        self.assertEqual(self.preview(form), "python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz")
        form._vars["n_directions"].set("500")
        form._vars["include_axis_views"].set(True)
        form._texts["direction"].insert("1.0", "-0.27,0.86,0.42")
        self.root.update()
        self.assertEqual(
            self.preview(form),
            "python analyze_projection_crossings.py examples/sc_W3_ALocal.xyz "
            "--n-directions 500 --include-axis-views --direction=-0.27,0.86,0.42",
        )

    def test_run_rejects_a_missing_input_file(self) -> None:
        form = self.make_form(analyze, ["examples/no_such_file.xyz", "--gui"])
        with mock.patch.object(form._messagebox, "showerror") as showerror:
            form.run()
        showerror.assert_called_once()
        self.assertIn("no_such_file.xyz", showerror.call_args[0][1])
        self.assertFalse(form.running)

    def test_run_executes_the_script_and_streams_its_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = os.path.join(tmp, "scan.csv")
            summary_path = os.path.join(tmp, "summary.txt")
            form = self.make_form(analyze, [
                "examples/sc_W3_ALocal.xyz", "--n-directions", "30", "--analysis-points", "200",
                "--progress-every", "10", "--csv-out", csv_path, "--summary-out", summary_path,
            ])
            with mock.patch.object(form._messagebox, "showerror", side_effect=AssertionError):
                form.run()
            self.assertTrue(form.running)
            deadline = time.monotonic() + 120
            while form.running and time.monotonic() < deadline:
                self.root.update()
                time.sleep(0.02)
            log = form.log_text.get("1.0", "end")
            self.assertEqual(form.returncode, 0, log)
            self.assertIn("analyzed 30/30 directions", log)
            self.assertIn("Finished successfully.", log)
            self.assertTrue(os.path.isfile(csv_path))
            with open(summary_path, encoding="utf-8") as handle:
                self.assertIn("Spherical scan directions: 30", handle.read())

    def test_stop_ends_a_running_script(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            form = self.make_form(analyze, [
                "examples/sc_W3_ALocal.xyz", "--n-directions", "200000", "--analysis-points", "0",
                "--engine", "numpy", "--csv-out", os.path.join(tmp, "scan.csv"),
                "--summary-out", os.path.join(tmp, "summary.txt"),
            ])
            form.run()
            self.assertTrue(form.running)
            form.stop()
            deadline = time.monotonic() + 30
            while form.running and time.monotonic() < deadline:
                self.root.update()
                time.sleep(0.02)
            self.assertFalse(form.running)
            self.assertNotEqual(form.returncode, 0)
            self.assertIn("Stopped.", form.log_text.get("1.0", "end"))


if __name__ == "__main__":
    unittest.main()
