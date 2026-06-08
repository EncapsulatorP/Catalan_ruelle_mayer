from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import mpmath as mp
import numpy as np

mp.mp.dps = 80

CATALAN = float(mp.catalan)
GOLDEN_RATIO = (1.0 + math.sqrt(5.0)) / 2.0
GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))


@dataclass(frozen=True)
class OutputPaths:
    root: Path
    figures: Path
    data: Path
    html: Path


def resolve_output_paths(root: Path, layout: str) -> OutputPaths:
    if layout not in {"canonical", "legacy"}:
        raise ValueError(f"Unsupported layout: {layout}")

    if layout == "legacy":
        paths = OutputPaths(root=root, figures=root, data=root, html=root / "html")
    else:
        paths = OutputPaths(
            root=root,
            figures=root / "figures",
            data=root / "data",
            html=root / "html",
        )

    for path in (paths.root, paths.figures, paths.data, paths.html):
        path.mkdir(parents=True, exist_ok=True)
    return paths


def li_s_i(s: float) -> complex:
    return mp.polylog(s, 1j)


def mg_channels(s: float) -> tuple[float, float]:
    value = li_s_i(s)
    denom = mp.pi ** s
    return float(mp.re(value) / denom), float(mp.im(value) / denom)


def signed_log_height(value: float, base: float = 10.0, eps: float = 1e-300) -> float:
    if abs(value) < eps:
        return 0.0
    return math.copysign(-math.log(abs(value), base), value)


def layer_type(s: float, channel: int) -> str:
    s_int = int(round(s))
    if channel == 0 and s_int % 2 == 0:
        return "rational"
    if channel == 1 and s_int % 2 == 1:
        return "rational"
    return "deep"


def layer_color(s: float, channel: int) -> str:
    return "#3498db" if layer_type(s, channel) == "rational" else "#d7263d"


def build_discrete_lattice(s_min: int = 2, s_max: int = 12) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for s in range(s_min, s_max + 1):
        a_value, b_value = mg_channels(s)
        for channel, value, name in (
            (0, a_value, "Re Li_s(i) / pi^s"),
            (1, b_value, "Im Li_s(i) / pi^s"),
        ):
            rho = -math.log10(abs(value)) if value else math.inf
            rows.append(
                {
                    "s": s,
                    "channel": channel,
                    "channel_name": name,
                    "value": value,
                    "rho": rho,
                    "signed_log_height": signed_log_height(value),
                    "type": layer_type(s, channel),
                    "color": layer_color(s, channel),
                }
            )
    return rows


def build_continuous_lattice(
    s_min: float = 2.0, s_max: float = 12.0, samples: int = 320
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for s in np.linspace(s_min, s_max, samples):
        a_value, b_value = mg_channels(float(s))
        for channel, value, name in (
            (0, a_value, "Re Li_s(i) / pi^s"),
            (1, b_value, "Im Li_s(i) / pi^s"),
        ):
            rho = -math.log10(abs(value)) if value else math.inf
            rows.append(
                {
                    "s": float(s),
                    "channel": channel,
                    "channel_name": name,
                    "value": value,
                    "rho": rho,
                    "signed_log_height": signed_log_height(value),
                }
            )
    return rows


def export_csv(rows: list[dict[str, object]], path: Path) -> None:
    if not rows:
        raise ValueError("Cannot export an empty table")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def rational_label(value: float, max_denominator: int = 2_000_000) -> str | None:
    fraction = Fraction(value).limit_denominator(max_denominator)
    if abs(float(fraction) - value) > 1e-12:
        return None
    sign = "-" if fraction.numerator < 0 else ""
    numerator = abs(fraction.numerator)
    return f"{sign}{numerator}/{fraction.denominator}"


def nested_radical(depth: int) -> float:
    value = 1.0
    for k in range(depth + 2, 2, -1):
        value = math.sqrt(1.0 + k * value)
    return value


def catalan_partial(depth: int) -> float:
    return (4.0 * CATALAN) / nested_radical(depth)


def stereographic_to_sphere(z: complex) -> tuple[float, float, float]:
    mod_sq = (z.real * z.real) + (z.imag * z.imag)
    denom = 1.0 + mod_sq
    x = 2.0 * z.real / denom
    y = 2.0 * z.imag / denom
    sphere_z = (mod_sq - 1.0) / denom
    return x, y, sphere_z


def pressure_proxy(
    t_values: np.ndarray, u_values: np.ndarray, k_max: int = 6000
) -> np.ndarray:
    k = np.arange(k_max, dtype=float)
    odd = 2.0 * k + 1.0
    exponent = u_values[..., None] - (2.0 * t_values[..., None])
    term_power = np.power(odd, exponent)
    return np.log(term_power.sum(axis=-1))


def first_n_primes(count: int) -> np.ndarray:
    if count < 1:
        return np.array([], dtype=int)

    if count < 6:
        bound = 15
    else:
        n = float(count)
        bound = int(n * (math.log(n) + math.log(math.log(n)))) + 10

    while True:
        sieve = np.ones(bound + 1, dtype=bool)
        sieve[:2] = False
        limit = int(math.sqrt(bound)) + 1
        for p in range(2, limit):
            if sieve[p]:
                sieve[p * p : bound + 1 : p] = False
        primes = np.flatnonzero(sieve)
        if len(primes) >= count:
            return primes[:count]
        bound *= 2


def primes_up_to(limit: int) -> np.ndarray:
    if limit < 2:
        return np.array([], dtype=int)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    boundary = int(math.sqrt(limit)) + 1
    for p in range(2, boundary):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = False
    return np.flatnonzero(sieve)


def half_integer_pi_distance(values: np.ndarray) -> np.ndarray:
    scaled = values / math.pi
    nearest = np.floor(scaled) + 0.5
    return np.abs(scaled - nearest)


def ti2(x: float) -> float:
    return float(mp.im(mp.polylog(2, 1j * x)))


def write_manifest(rows: Iterable[dict[str, object]], path: Path) -> None:
    rows = list(rows)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
