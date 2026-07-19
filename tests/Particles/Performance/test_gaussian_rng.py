#!/usr/bin/env python3
"""Statistical checks for shared polar truncated-Gaussian sampling."""

from __future__ import annotations

import math
import random
import unittest


def polar_normal_pair(u: float, v: float) -> tuple[float, float] | None:
    x = 2.0 * u - 1.0
    y = 2.0 * v - 1.0
    s = x * x + y * y
    if not (0.0 < s < 1.0):
        return None
    factor = math.sqrt(-2.0 * math.log(s) / s)
    return x * factor, y * factor


def scale_truncated(z: float, mu: float, sigma: float, trunc_limit: float, trunc_factor: float) -> float | None:
    bound = trunc_limit / trunc_factor
    if abs(z) <= bound:
        return mu + sigma * trunc_factor * z
    return None


class PolarSampler:
    def __init__(self, trunc_limit: float = 2.0, trunc_factor: float = 1.05) -> None:
        self.trunc_limit = trunc_limit
        self.trunc_factor = trunc_factor
        self.spare_ready = False
        self.spare_std = 0.0
        self._rng = random.Random()

    def seed(self, value: int) -> None:
        self._rng.seed(value)
        self.spare_ready = False
        self.spare_std = 0.0

    def sample(self, mu: float = 0.0, sigma: float = 1.0) -> float:
        while True:
            if self.spare_ready:
                z = self.spare_std
                self.spare_ready = False
            else:
                while True:
                    pair = polar_normal_pair(self._rng.random(), self._rng.random())
                    if pair is not None:
                        z, spare = pair
                        self.spare_std = spare
                        self.spare_ready = True
                        break
            value = scale_truncated(z, mu, sigma, self.trunc_limit, self.trunc_factor)
            if value is not None:
                return value


class GaussianRngTests(unittest.TestCase):
    def test_samples_are_finite_and_truncated(self) -> None:
        sampler = PolarSampler()
        sampler.seed(104729)
        samples = [sampler.sample(sigma=0.25) for _ in range(20_000)]
        self.assertTrue(all(math.isfinite(v) for v in samples))
        limit = 2.0 * 0.25
        self.assertTrue(all(abs(v) <= limit + 1e-12 for v in samples))

    def test_mean_and_variance_with_wide_truncation(self) -> None:
        sampler = PolarSampler(trunc_limit=4.0, trunc_factor=1.0)
        sampler.seed(130363)
        sigma = 0.5
        samples = [sampler.sample(sigma=sigma) for _ in range(80_000)]
        mean = sum(samples) / len(samples)
        variance = sum((v - mean) ** 2 for v in samples) / (len(samples) - 1)
        self.assertLess(abs(mean), 0.02)
        self.assertLess(abs(math.sqrt(variance) - sigma), 0.03)

    def test_fixed_seed_reproducible(self) -> None:
        a = PolarSampler()
        b = PolarSampler()
        a.seed(155921)
        b.seed(155921)
        self.assertEqual([a.sample() for _ in range(256)], [b.sample() for _ in range(256)])

    def test_pure_transforms_from_fixed_uniforms(self) -> None:
        # u,v chosen so (2u-1)^2+(2v-1)^2 is strictly inside the unit disk.
        pair = polar_normal_pair(0.3, 0.4)
        self.assertIsNotNone(pair)
        assert pair is not None
        scaled = scale_truncated(pair[0], 0.0, 1.0, 3.0, 1.0)
        self.assertIsNotNone(scaled)
        self.assertTrue(math.isfinite(scaled))


if __name__ == "__main__":
    unittest.main()
