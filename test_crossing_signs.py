#!/usr/bin/env python3
"""Regression tests for crossing counts and the Gauss-consistent sign."""

import os
import unittest

import numpy as np

from analyze_projection_crossings import (
    count_projection_crossings,
    fibonacci_sphere,
    read_xyz,
    resample_closed_polyline,
)


HERE = os.path.dirname(os.path.abspath(__file__))
XYZ = os.path.join(HERE, "examples", "sc_W3_ALocal.xyz")


class CrossingSignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.full = read_xyz(XYZ)

    def test_coordinate_views_have_positive_sign(self) -> None:
        self.assertEqual(count_projection_crossings(self.full, [0, -1, 0])[:2], (3, 3))
        self.assertEqual(count_projection_crossings(self.full, [1, 0, 0])[:2], (4, 4))

    def test_viewing_from_opposite_side_preserves_directional_writhe(self) -> None:
        normal = np.array([0.2494, 0.8564, 0.4521])
        self.assertEqual(
            count_projection_crossings(self.full, normal)[:2],
            count_projection_crossings(self.full, -normal)[:2],
        )

    def test_full_cancellation_view(self) -> None:
        self.assertEqual(
            count_projection_crossings(self.full, [0.1941, 0.2748, 0.9417])[:2],
            (6, 0),
        )

    def test_spherical_mean_is_positive(self) -> None:
        points = resample_closed_polyline(self.full, 500)
        signed = [count_projection_crossings(points, n)[1] for n in fibonacci_sphere(200)]
        self.assertGreater(float(np.mean(signed)), 2.8)
        self.assertLess(float(np.mean(signed)), 3.2)


if __name__ == "__main__":
    unittest.main()
