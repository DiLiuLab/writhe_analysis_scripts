#!/usr/bin/env python3
"""Regression tests for the data-driven parts of compare_projection_scans.py."""

import os
import tempfile
import unittest

import compare_projection_scans as compare


HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "examples")
COMPARISONS = os.path.join(EXAMPLES, "conformation_comparison_figures")


def matched(scan_a: str, scan_b: str):
    return compare.align_rows(
        compare.read_fibonacci_rows(os.path.join(EXAMPLES, scan_a)),
        compare.read_fibonacci_rows(os.path.join(EXAMPLES, scan_b)),
    )


class CompareSummaryTests(unittest.TestCase):
    def test_positive_pair_keeps_the_original_rows_and_title(self) -> None:
        data = matched("sc_W3_ALocal_projection_scan.csv", "sc_W3_R11_projection_scan.csv")
        self.assertEqual(compare.highlighted_values(data), ([3], [3, 4]))
        self.assertEqual(
            compare.coordinate_view_title(data),
            "Coordinate-plane views of two Wr approximately +3 conformations",
        )

    def test_positive_pair_summary_matches_the_committed_file(self) -> None:
        data = matched("sc_W3_ALocal_projection_scan.csv", "sc_W3_R11_projection_scan.csv")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "summary.txt")
            compare.write_summary(path, data, "ALocal", "R11")
            with open(path, encoding="utf-8") as new, open(
                os.path.join(COMPARISONS, "ALocal_vs_R11_comparison_summary.txt"), encoding="utf-8"
            ) as committed:
                self.assertEqual(new.read(), committed.read())

    def test_negative_pair_uses_its_own_values(self) -> None:
        data = matched("sc_L714_Wm2_Logo_projection_scan.csv", "sc_L714_Wm2_R13_ABend_projection_scan.csv")
        self.assertEqual(compare.highlighted_values(data), ([2], [-2, -3]))
        self.assertEqual(
            compare.coordinate_view_title(data),
            "Coordinate-plane views of two Wr approximately -2 conformations",
        )

    def test_different_writhes_get_a_neutral_title(self) -> None:
        data = matched("sc_W3_ALocal_projection_scan.csv", "sc_L714_Wm2_Logo_projection_scan.csv")
        self.assertEqual(compare.coordinate_view_title(data), "Coordinate-plane views of two conformations")


if __name__ == "__main__":
    unittest.main()
