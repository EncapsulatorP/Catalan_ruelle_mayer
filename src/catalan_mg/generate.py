from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from .core import (
    CATALAN,
    GOLDEN_ANGLE,
    GOLDEN_RATIO,
    build_continuous_lattice,
    build_discrete_lattice,
    catalan_partial,
    export_csv,
    first_n_primes,
    half_integer_pi_distance,
    layer_type,
    mg_channels,
    nested_radical,
    pressure_proxy,
    primes_up_to,
    rational_label,
    resolve_output_paths,
    signed_log_height,
    stereographic_to_sphere,
    ti2,
    write_manifest,
)

warnings.filterwarnings(
    "ignore",
    message="constrained_layout not applied because axes sizes collapsed to zero",
)


def _save_figure(fig: plt.Figure, path: Path, dpi: int = 220) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _record(records: list[dict[str, object]], category: str, path: Path) -> None:
    records.append(
        {
            "category": category,
            "relative_path": path.as_posix(),
        }
    )


def _plot_checkerboard(rows: list[dict[str, object]], path: Path) -> None:
    fig = plt.figure(figsize=(11.5, 8.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    for row in rows:
        ax.scatter(
            row["s"],
            row["channel"],
            row["signed_log_height"],
            s=95,
            color=row["color"],
            edgecolor="black",
            linewidth=0.5,
            alpha=0.95,
        )
        ax.text(
            row["s"],
            row["channel"],
            row["signed_log_height"],
            f"s={row['s']}",
            fontsize=8,
            color=row["color"],
        )

    for channel in (0, 1):
        sheet = [row for row in rows if row["channel"] == channel]
        ax.plot(
            [row["s"] for row in sheet],
            [row["channel"] for row in sheet],
            [row["signed_log_height"] for row in sheet],
            alpha=0.35,
        )

    ax.set_title(r"$\mathcal{M}_G$ Two-Sheet Signed-Log Lattice", fontsize=16)
    ax.set_xlabel("weight s")
    ax.set_ylabel("channel: 0=Re, 1=Im")
    ax.set_zlabel(r"sign(value) * -log10 |value|")
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#3498db",
                markeredgecolor="black",
                label="known rational-over-pi layer",
                markersize=10,
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#d7263d",
                markeredgecolor="black",
                label="deep / Catalan-like layer",
                markersize=10,
            ),
        ],
        loc="upper left",
    )

    _save_figure(fig, path)


def _plot_double_helix(rows: list[dict[str, object]], path: Path) -> None:
    fig = plt.figure(figsize=(11.5, 8.8), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    for channel, color, label in (
        (0, "#3498db", "Re / zeta channel"),
        (1, "#d7263d", "Im / beta-Catalan channel"),
    ):
        sheet = [row for row in rows if row["channel"] == channel]
        xs: list[float] = []
        ys: list[float] = []
        zs: list[float] = []

        for row in sheet:
            theta = (math.pi / 2.0) * row["s"] + (math.pi * channel)
            rho = row["rho"]
            xs.append(rho * math.cos(theta))
            ys.append(rho * math.sin(theta))
            zs.append(row["s"])

        ax.plot(xs, ys, zs, color=color, linewidth=2.0, label=label)
        ax.scatter(xs, ys, zs, color=color, s=8, alpha=0.55)

    ax.set_title(r"$\mathcal{M}_G$ Double-Helix Projection", fontsize=16)
    ax.set_xlabel(r"$\rho \cos \theta$")
    ax.set_ylabel(r"$\rho \sin \theta$")
    ax.set_zlabel("weight s")
    ax.legend()

    _save_figure(fig, path)


def _plot_prime_style_funnel(
    rows: list[dict[str, object]], path: Path, title: str
) -> None:
    fig = plt.figure(figsize=(12.0, 9.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    x_values: list[float] = []
    y_values: list[float] = []
    z_values: list[float] = []
    colors: list[str] = []

    for row in rows:
        theta = (math.pi / 2.0) * row["s"] + (math.pi * row["channel"])
        rho = row["rho"]
        x_values.append(rho * math.sin(theta) * math.cos(theta))
        y_values.append(rho * (math.sin(theta) ** 2))
        z_values.append(row["s"] * math.cos(theta))
        colors.append("#3498db" if row["channel"] == 0 else "#d7263d")

    ax.scatter(
        x_values,
        y_values,
        z_values,
        c=colors,
        s=12,
        alpha=0.75,
        edgecolor="none",
    )
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(r"$x = \rho \sin \theta \cos \theta$")
    ax.set_ylabel(r"$y = \rho \sin^2 \theta$")
    ax.set_zlabel(r"$z = s \cos \theta$")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#3498db", label="Re / zeta channel", markersize=9),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#d7263d", label="Im / beta-Catalan channel", markersize=9),
        ],
        loc="upper left",
    )

    _save_figure(fig, path)


def _plot_catalan_apery_heatmap(rows: list[dict[str, object]], path: Path) -> None:
    fig = plt.figure(figsize=(12.0, 8.8), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    type_values = np.array([0 if row["type"] == "rational" else 1 for row in rows])
    color_map = ListedColormap(["#4292e8", "#e53459"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], color_map.N)

    scatter = ax.scatter(
        [row["s"] for row in rows],
        [row["channel"] for row in rows],
        [math.log10(abs(row["value"])) for row in rows],
        c=type_values,
        cmap=color_map,
        norm=norm,
        s=85,
        edgecolor="#303030",
        linewidth=0.6,
        alpha=0.95,
    )
    for channel in (0, 1):
        sheet = [row for row in rows if row["channel"] == channel]
        ax.plot(
            [row["s"] for row in sheet],
            [row["channel"] for row in sheet],
            [math.log10(abs(row["value"])) for row in sheet],
            color="#8f8f8f",
            alpha=0.35,
        )

    ax.set_title(
        r"$M_G(s)$ -- Catalan-Apery Checkerboard"
        "\n(Blue = Rational layer, Crimson = Deep layer)",
        fontsize=15,
    )
    ax.set_xlabel("Weight s (integer >= 2)")
    ax.set_ylabel("Channel")
    ax.set_zlabel(r"$\log_{10} |$channel value$|$")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Re / pi^s", "Im / pi^s"])
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.12)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Rational", "Deep"])
    cbar.set_label("Layer Type")

    _save_figure(fig, path)


def _plot_dark_mg_matrix(rows: list[dict[str, object]], path: Path) -> None:
    fig = plt.figure(figsize=(10.0, 8.5), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("#0f131a")
    ax.set_facecolor("#0f131a")

    x_values = []
    y_values = []
    z_values = []
    colors = []
    for row in rows:
        if row["channel"] == 0:
            x_values.append(row["s"])
            y_values.append(row["value"])
            z_values.append(0.0)
            colors.append("#5ec7ff" if row["type"] == "rational" else "#ff7d50")
        else:
            x_values.append(row["s"])
            y_values.append(0.0)
            z_values.append(row["value"])
            colors.append("#5ec7ff" if row["type"] == "rational" else "#ff7d50")

    ax.scatter(x_values, y_values, z_values, c=colors, s=170, edgecolor="#e8ebef", linewidth=1.0)

    re_rows = [row for row in rows if row["channel"] == 0]
    im_rows = [row for row in rows if row["channel"] == 1]
    ax.plot(
        [row["s"] for row in re_rows],
        [row["value"] for row in re_rows],
        [0.0 for _ in re_rows],
        color="#2f8cff",
        alpha=0.45,
        linewidth=2.0,
    )
    ax.plot(
        [row["s"] for row in im_rows],
        [0.0 for _ in im_rows],
        [row["value"] for row in im_rows],
        color="#ff6b3d",
        alpha=0.45,
        linewidth=2.0,
    )

    for row in rows:
        if row["type"] == "rational":
            label = rational_label(row["value"])
            if not label:
                continue
            if row["channel"] == 0:
                ax.text(row["s"], row["value"], 0.004, label, color="#adff76", fontsize=10)
            else:
                ax.text(row["s"], 0.001, row["value"], label, color="#adff76", fontsize=10)
        else:
            if row["channel"] == 0:
                ax.text(row["s"], row["value"], 0.001, f"s={row['s']}", color="#5ec7ff", fontsize=10)
            else:
                ax.text(row["s"], 0.001, row["value"], f"s={row['s']}", color="#ff7d50", fontsize=10)

    ax.set_title(
        "Mirror-Catalan Matrix\n"
        r"$\mathcal{M}_G(s) = (\mathrm{Re}\,\mathrm{Li}_s(i)/\pi^s,\ \mathrm{Im}\,\mathrm{Li}_s(i)/\pi^s)$",
        color="#f2f5f8",
        pad=20,
    )
    ax.set_xlabel("s  (weight / level)", color="#f2f5f8")
    ax.set_ylabel(r"Re Li_s(i) / pi^s", color="#5ec7ff")
    ax.set_zlabel(r"Im Li_s(i) / pi^s", color="#ff7d50")
    ax.tick_params(colors="#f2f5f8")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((0.08, 0.10, 0.13, 1.0))
    ax.legend(
        handles=[
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#5ec7ff", label="rational channel", markersize=12),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#ff7d50", label="deep channel", markersize=12),
        ],
        facecolor="#1b2230",
        edgecolor="#707788",
        labelcolor="#f2f5f8",
        loc="upper right",
    )

    _save_figure(fig, path)


def _plot_prime_projection(path: Path) -> None:
    primes = first_n_primes(3000).astype(float)
    x_values = primes * np.sin(primes) * np.cos(primes)
    y_values = primes * np.sin(primes) ** 2
    z_values = primes * np.cos(primes)

    fig = plt.figure(figsize=(12.0, 7.6), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    scatter = ax.scatter(
        x_values,
        y_values,
        z_values,
        c=primes,
        cmap="plasma",
        s=5,
        alpha=0.85,
    )
    ax.set_title(r"$x = p \sin(p)\cos(p),\ y = p\sin^2(p),\ z = p\cos(p)$")
    ax.set_xlabel(r"$x = p \sin(p)\cos(p)$")
    ax.set_ylabel(r"$y = p\sin^2(p)$")
    ax.set_zlabel(r"$z = p\cos(p)$")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.1)
    cbar.set_label("Prime Value (p)")

    _save_figure(fig, path)


def _stable_log_abs_tanh_residual(x_values: np.ndarray) -> np.ndarray:
    abs_x = np.abs(x_values)
    direct = np.log(np.maximum(1.0 - np.abs(np.tanh(x_values)), 1e-300))
    asymptotic = math.log(2.0) - (2.0 * abs_x)
    return np.where(abs_x > 18.0, asymptotic, direct)


def _plot_tanh_diophantine(path: Path) -> None:
    primes = primes_up_to(200000).astype(float)
    distance = half_integer_pi_distance(primes)
    x_values = primes * np.cos(primes)
    log_residual = _stable_log_abs_tanh_residual(x_values)
    diophantine_cost = 2.0 * math.pi * primes * distance
    reference = diophantine_cost - math.log(2.0)

    anomalies = np.argsort(log_residual)[-15:]
    anomalies = anomalies[np.argsort(primes[anomalies])]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)

    mask = distance < 0.05
    axes[0].scatter(
        diophantine_cost[mask],
        -log_residual[mask],
        s=4,
        alpha=0.45,
        color="#2d72b8",
    )
    x_line = np.linspace(0.0, float(diophantine_cost[mask].max()) * 1.02, 300)
    axes[0].plot(x_line, x_line - math.log(2.0), color="#a91e1e", linewidth=1.2)
    axes[0].set_title(
        "Saturation residual vs half-integer-pi distance\n"
        "(primes with d(p) < 0.05)"
    )
    axes[0].set_xlabel(r"$2\pi p \cdot d(p)$  (Diophantine cost)")
    axes[0].set_ylabel(r"$-\log(1 - |\tanh(p\cos p)|)$")
    axes[0].grid(alpha=0.2)

    axes[1].scatter(primes, log_residual, s=2, alpha=0.18, color="#676767")
    axes[1].scatter(primes[anomalies], log_residual[anomalies], s=45, color="#e84b4b", label="top 15 anomalies")
    for index in anomalies[:8]:
        axes[1].text(primes[index], log_residual[index] + 2000.0, f"{int(primes[index])}", color="#9e2828", fontsize=8)
    axes[1].set_title(
        "Anomalous primes:\n"
        r"those where $\tanh(p\cos p)$ fails to saturate"
    )
    axes[1].set_xlabel("prime p")
    axes[1].set_ylabel(r"$\log r(p)$")
    axes[1].grid(alpha=0.2)
    axes[1].legend(loc="lower right", frameon=False)

    _save_figure(fig, path)


def _plot_ti2_primes(path: Path) -> None:
    primes = primes_up_to(500)
    ti2_values = np.array([ti2(float(p)) for p in primes], dtype=float)
    left_values = ti2_values - CATALAN
    log_asymptotic = (math.pi / 2.0) * np.log(primes)
    regularized = np.array([ti2(1.0 / float(p)) for p in primes], dtype=float)
    leading_term = 1.0 / primes
    corrected_term = leading_term - (1.0 / (9.0 * (primes ** 3)))

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.5), constrained_layout=True)

    axes[0].scatter(primes, left_values, s=18, color="#1f66ad", label=r"$\mathrm{Ti}_2(p) - G$")
    axes[0].plot(primes, log_asymptotic, linestyle="--", color="#b63232", label=r"$(\pi/2)\log p$")
    axes[0].plot(primes, log_asymptotic - CATALAN, linestyle=":", color="#17a276", label=r"$(\pi/2)\log p - G$")
    axes[0].set_title(r"$\mathrm{Ti}_2(p) - G$: diverges like $(\pi/2)\log p$")
    axes[0].set_xlabel("prime p")
    axes[0].set_ylabel("value")
    axes[0].grid(alpha=0.22)
    axes[0].legend(frameon=False)

    axes[1].scatter(primes, regularized, s=18, color="#1f66ad", label=r"$\mathrm{Ti}_2(1/p)$")
    axes[1].plot(primes, leading_term, linestyle="--", color="#b63232", label=r"$1/p$  (leading term)")
    axes[1].plot(primes, corrected_term, linestyle=":", color="#17a276", label=r"$1/p - 1/(9p^3)$")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_title(
        r"$\mathrm{Ti}_2(p) - (\pi/2)\log p - G = \mathrm{Ti}_2(1/p)$"
        "\nso the regularized object decays like 1/p"
    )
    axes[1].set_xlabel("prime p")
    axes[1].set_ylabel(r"$\mathrm{Ti}_2(1/p)$")
    axes[1].grid(alpha=0.22, which="both")
    axes[1].legend(frameon=False)

    _save_figure(fig, path)


def _plot_pressure_slice(path: Path) -> None:
    u_values = np.linspace(-0.5, 1.3, 260)
    t_value = np.full_like(u_values, 1.0)
    pressure = pressure_proxy(t_value, u_values, k_max=5000)

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.plot(u_values, pressure, color="#1b40f0", linewidth=2.0)
    ax.axvline(1.0, color="red", linestyle="--", label="Theoretical boundary u=2t-1 = 1.0")
    ax.set_title("Pressure at t=1.0")
    ax.set_xlabel("u")
    ax.set_ylabel("P(1.0, u)")
    ax.grid(alpha=0.25)
    ax.legend()

    _save_figure(fig, path)


def _plot_pressure_phase_diagram(path: Path) -> None:
    t_values = np.linspace(0.55, 1.4, 120)
    u_values = np.linspace(-0.5, 1.4, 120)
    t_grid, u_grid = np.meshgrid(t_values, u_values)
    pressure = pressure_proxy(t_grid, u_grid, k_max=2500)

    finite_cutoff = []
    for t_value in t_values:
        u_line = np.linspace(-0.5, 1.4, 240)
        large = pressure_proxy(np.full_like(u_line, t_value), u_line, k_max=5000)
        small = pressure_proxy(np.full_like(u_line, t_value), u_line, k_max=500)
        gap = np.abs(large - small)
        crossing = np.where(gap > 0.18)[0]
        finite_cutoff.append(float(u_line[crossing[0]]) if crossing.size else math.nan)
    finite_cutoff = np.array(finite_cutoff)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    contour = axes[0].contourf(t_values, u_values, pressure, levels=18, cmap="viridis")
    fig.colorbar(contour, ax=axes[0], fraction=0.05, pad=0.04)
    axes[0].set_title(r"Pressure $P(-t\log|T'| + u\log(2k+1))$")
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("u")

    axes[1].plot(t_values, 2.0 * t_values - 1.0, color="red", label="Theoretical: u = 2t-1")
    axes[1].plot(t_values, finite_cutoff, color="blue", linestyle="--", label="Finite-cutoff proxy")
    axes[1].set_title("Phase boundary")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("u")
    axes[1].legend(loc="upper left")
    axes[1].grid(alpha=0.25)

    _save_figure(fig, path)


def _build_zded_lattice(max_k: int = 12) -> list[dict[str, object]]:
    from mpmath import mpf, zeta

    lattice: list[dict[str, object]] = []
    for k in range(max_k):
        s_value = mpf(k + 2)
        layer = "rational" if k % 2 == 0 else "irrational"
        lattice.append(
            {
                "k": k,
                "s": s_value,
                "layer": layer,
                "exact": zeta(s_value),
            }
        )
    return lattice


def _zded_thermo_at_s(
    s_value: float,
    n_start: int = 100,
    n_steps: int = 80,
    beta0: float = 1.0,
    mu0: float = 0.5,
    seed: int = 42,
) -> tuple[float, float, float, float]:
    from mpmath import log, mpf

    rng = np.random.default_rng(seed + int(round(float(s_value) * 100)))
    word_k = rng.geometric(p=0.5, size=n_steps)

    rho = mpf(1.0)
    frac = mpf(0.0)
    action = mpf(0.0)
    n_t = mpf(n_start)
    mu = mu0 * s_value
    beta = beta0 / (s_value / 2.0) if s_value > 0 else beta0

    for word in word_k:
        r_value = mpf(3.0)
        b_value = mpf(1 << int(word))
        n_next = (r_value * n_t + 1.0) / b_value
        rho_new = (r_value / b_value) * rho
        frac_new = (r_value / b_value) * frac + (1.0 / b_value)

        if n_t > 0 and n_next > 0:
            action += beta * word + mu * log(n_t / n_next)
        else:
            action += beta * word

        rho = rho_new
        frac = frac_new
        n_t = n_next

    one_minus_rho = max(1.0 - rho, mpf("1e-30"))
    reduced = frac % one_minus_rho
    lock_error = abs(reduced) / one_minus_rho
    return float(rho), float(action), float(lock_error), float(rho)


def _model_value(rho: float, action: float, lock_error: float, r_over_b: float) -> float:
    return (1.0 / (1.0 + lock_error)) * math.exp(-0.01 * action) * (1.0 + 0.5 * r_over_b)


def _calibrate_zded(
    lattice: list[dict[str, object]],
    n_trials: int = 120,
    seed: int = 7,
) -> tuple[tuple[float, float], float]:
    rng = np.random.default_rng(seed)
    rational_points = [point for point in lattice if point["layer"] == "rational"]
    best_error = float("inf")
    best_params = (0.5, 1.0)

    for _ in range(n_trials):
        mu0 = rng.uniform(0.1, 2.0)
        beta0 = rng.uniform(0.5, 3.0)
        total_error = 0.0
        for point in rational_points:
            rho, action, lock_error, rb = _zded_thermo_at_s(float(point["s"]), mu0=mu0, beta0=beta0)
            total_error += abs(_model_value(rho, action, lock_error, rb) - float(point["exact"]))
        average_error = total_error / len(rational_points)
        if average_error < best_error:
            best_error = average_error
            best_params = (mu0, beta0)

    return best_params, best_error


def _probe_zded(
    lattice: list[dict[str, object]], mu0: float, beta0: float
) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    for point in lattice:
        if point["layer"] != "irrational":
            continue
        rho, action, lock_error, rb = _zded_thermo_at_s(float(point["s"]), mu0=mu0, beta0=beta0)
        hypothesis = _model_value(rho, action, lock_error, rb)
        true_value = float(point["exact"])
        rel_error = abs(hypothesis - true_value) / true_value if true_value else 0.0
        results.append(
            {
                "s": float(point["s"]),
                "v_hypothesis": hypothesis,
                "true_zeta": true_value,
                "rel_error": rel_error,
                "lock_error": lock_error,
                "L_action": action,
                "rho": rho,
            }
        )
    return results


def _plot_zded_calibration(fig_path: Path, csv_path: Path) -> None:
    lattice = _build_zded_lattice(max_k=12)
    (mu0_best, beta0_best), calibration_error = _calibrate_zded(lattice)
    irrational_results = _probe_zded(lattice, mu0_best, beta0_best)

    export_csv(irrational_results, csv_path)

    calibration_label = (
        f"Toy calibration mu0={mu0_best:.3f}, beta0={beta0_best:.3f}, residual={calibration_error:.2e}"
    )

    fig, ax = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    rational_s = [float(point["s"]) for point in lattice if point["layer"] == "rational"]
    rational_zeta = [float(point["exact"]) for point in lattice if point["layer"] == "rational"]
    irrational_s = [result["s"] for result in irrational_results]
    irrational_hyp = [result["v_hypothesis"] for result in irrational_results]
    irrational_true = [result["true_zeta"] for result in irrational_results]

    ax.plot(rational_s, rational_zeta, "o-", color="#2a6fdb", label="Rational layers (exact Q(pi))", markersize=7)
    ax.plot(irrational_s, irrational_hyp, "s--", color="#d1495b", label="Irrational hypothesis (model)", markersize=7)
    ax.plot(irrational_s, irrational_true, "^", color="#2f9e44", label="True zeta reference", markersize=6, alpha=0.8)
    ax.set_xlabel("s (Mg lattice layer)")
    ax.set_ylabel("Effective value / zeta(s)")
    ax.set_title(
        "Exploratory ZDED Mg lattice toy calibration"
        f"\n{calibration_label}"
    )
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend()

    _save_figure(fig, fig_path)


def _plot_catalan_radical_riemann(path: Path) -> None:
    depths = np.arange(1, 35)
    radicals = np.array([nested_radical(int(depth)) for depth in depths])
    g_partials = np.array([catalan_partial(int(depth)) for depth in depths])
    g_errors = np.abs(g_partials - CATALAN)

    fig = plt.figure(figsize=(15.0, 8.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.0, 1.0])

    ax_sphere = fig.add_subplot(grid[:, 0], projection="3d")
    ax_arc = fig.add_subplot(grid[0, 1])
    ax_conv = fig.add_subplot(grid[1, 1])

    fig.patch.set_facecolor("#0b0c14")
    ax_sphere.set_facecolor("#0f101a")
    for axis in (ax_arc, ax_conv):
        axis.set_facecolor("#0f101a")

    u = np.linspace(0.0, 2.0 * math.pi, 60)
    v = np.linspace(0.0, math.pi, 30)
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones_like(u), np.cos(v))
    ax_sphere.plot_wireframe(sphere_x, sphere_y, sphere_z, color="#262c44", linewidth=0.35, alpha=0.22)

    colors = plt.cm.plasma(np.linspace(0.15, 0.95, len(depths)))
    sphere_points = []
    for depth, radical in zip(depths, radicals):
        complex_point = (math.pi / radical) * complex(math.cos(depth * GOLDEN_ANGLE), math.sin(depth * GOLDEN_ANGLE))
        sphere_points.append(stereographic_to_sphere(complex_point))
    sphere_points = np.array(sphere_points)

    ax_sphere.plot(
        sphere_points[:, 0],
        sphere_points[:, 1],
        sphere_points[:, 2],
        color="#8338ec",
        linewidth=2.4,
        alpha=0.85,
    )
    ax_sphere.scatter(
        sphere_points[:, 0],
        sphere_points[:, 1],
        sphere_points[:, 2],
        c=colors,
        s=np.linspace(42, 12, len(depths)),
        edgecolor="#ffffff",
        linewidth=0.4,
        alpha=0.9,
    )
    ax_sphere.scatter([0], [0], [-1], marker="*", s=280, color="#f1c453", edgecolor="#d9dfe8", linewidth=1.3)
    for index in (0, 3, 9, 19, 33):
        ax_sphere.text(
            sphere_points[index, 0],
            sphere_points[index, 1],
            sphere_points[index, 2] + 0.08,
            f"N={depths[index]}",
            color=colors[index],
            fontsize=10,
        )
    ax_sphere.text2D(
        0.02,
        0.98,
        "Phyllotactic stereographic projection on the Riemann sphere\n"
        r"$z_N = (\pi / R_N)e^{iN\varphi}$,  radius scales with the radical truncation",
        transform=ax_sphere.transAxes,
        color="#ebedf3",
        fontsize=14,
        va="top",
    )
    for axis in (ax_sphere.xaxis, ax_sphere.yaxis, ax_sphere.zaxis):
        axis.pane.set_facecolor((0.06, 0.07, 0.1, 1.0))
        axis.pane.set_edgecolor((0.18, 0.20, 0.28, 0.55))
    ax_sphere.set_axis_off()
    ax_sphere.view_init(elev=23, azim=-64)

    ax_arc.text(
        0.5,
        0.98,
        r"$R_N = \sqrt{1 + 3\sqrt{1 + 4\sqrt{1 + 5\sqrt{\cdots}}}} \to 4$",
        transform=ax_arc.transAxes,
        ha="center",
        va="top",
        color="#ebedf3",
        fontsize=16,
    )
    for idx, label in enumerate(range(3, 10)):
        radius = 0.8 - (idx * 0.1)
        theta = np.linspace(0.15 * math.pi, 0.85 * math.pi, 200)
        arc_x = radius * np.cos(theta)
        arc_y = radius * np.sin(theta)
        ax_arc.plot(arc_x, arc_y, color=plt.cm.viridis(0.05 + idx / 9.0), linewidth=2.4)
        ax_arc.text(0.0, radius + 0.03, str(label), color=plt.cm.viridis(0.05 + idx / 9.0), ha="center", fontsize=11)
    ax_arc.text(-0.55, -0.25, f"depth N=20: R = {nested_radical(20):.7f}", color="#a6a9ff", fontsize=12, family="monospace")
    ax_arc.text(-0.55, -0.42, f"limit:      R = {4.0:.7f}", color="#f5d76e", fontsize=12, family="monospace")
    ax_arc.set_xlim(-1.0, 1.0)
    ax_arc.set_ylim(-0.5, 1.0)
    ax_arc.axis("off")

    ax_conv.semilogy(depths, g_errors, "o-", color="#f78c6b", markerfacecolor="#c93f87", markeredgecolor="#f2e6d8")
    ax_conv.semilogy(depths, 0.5 * np.power(0.25, depths - 1), "--", color="#c28d3d", label="0.5 * 4^-(N-1) reference")
    ax_conv.set_title("Geometric convergence to Catalan", color="#ebedf3")
    ax_conv.set_xlabel("truncation depth N", color="#ebedf3")
    ax_conv.set_ylabel(r"$|G_N - G|$", color="#ebedf3")
    ax_conv.grid(alpha=0.14)
    ax_conv.legend(facecolor="#121726", edgecolor="#616777", labelcolor="#ebedf3")
    ax_conv.tick_params(colors="#c8cce3")
    for spine in ax_conv.spines.values():
        spine.set_color("#41475c")

    fig.suptitle(
        "Convergence of the Ramanujan-style radical board toward Catalan",
        color="#ebedf3",
        fontsize=16,
        y=0.975,
    )

    _save_figure(fig, path)


def _plot_ramanujan_spiral_svg(path: Path) -> None:
    depths = np.arange(1, 31)
    radicals = np.array([nested_radical(int(depth)) for depth in depths])
    g_partials = np.array([catalan_partial(int(depth)) for depth in depths])
    radical_gap = np.abs(4.0 - radicals)
    radii = 240.0 * np.power(radical_gap / radical_gap.max(), 0.55)
    theta = depths * 0.55
    x_values = radii * np.cos(theta)
    y_values = radii * np.sin(theta)

    fig = plt.figure(figsize=(8.5, 9.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.1])
    ax_spiral = fig.add_subplot(grid[0, 0])
    ax_table = fig.add_subplot(grid[1, 0])

    fig.patch.set_facecolor("#f6f4ee")
    ax_spiral.set_facecolor("#0f101a")
    ax_table.set_facecolor("#f6f4ee")

    circle_radii = [60, 110, 160, 210, 250]
    for radius in circle_radii:
        ax_spiral.add_artist(plt.Circle((0, 0), radius, fill=False, color="#8b90a0", alpha=0.22, linewidth=0.8))

    ax_spiral.plot(x_values, y_values, color="#7f77dd", linewidth=1.8, alpha=0.92)
    color_values = plt.cm.Purples(np.linspace(0.85, 0.25, len(depths)))
    sizes = np.linspace(70, 15, len(depths))
    ax_spiral.scatter(x_values, y_values, s=sizes, c=color_values, edgecolor="#ffffff", linewidth=0.4)
    ax_spiral.scatter([0], [0], s=90, color="#1d9e75", edgecolor="#ffffff", linewidth=0.8, zorder=5)
    ax_spiral.text(7, 8, "G", color="#1d9e75", fontsize=11, weight="bold")
    for index in (0, 1, 2, 3, 6, 9, 14):
        ax_spiral.text(x_values[index] + 8, y_values[index] + 8, f"d={depths[index]}", color="#ece8ff", fontsize=9)

    ax_spiral.set_title("Ramanujan recursion -> Catalan", color="#f8f9fb", fontsize=16)
    ax_spiral.text(
        0.5,
        0.95,
        "Polar log-spiral: angle = depth, radius = distance from the fixed point 4",
        color="#c4c0b6",
        transform=ax_spiral.transAxes,
        ha="center",
        fontsize=10,
    )
    ax_spiral.set_xlim(-270, 270)
    ax_spiral.set_ylim(-270, 270)
    ax_spiral.set_aspect("equal")
    ax_spiral.axis("off")

    ax_table.axis("off")
    ax_table.text(0.02, 0.88, "depth d", fontsize=11, weight="bold")
    ax_table.text(0.16, 0.88, "nested radical", fontsize=11, weight="bold")
    ax_table.text(0.46, 0.88, "G partial = 4G / radical", fontsize=11, weight="bold")
    ax_table.text(0.78, 0.88, "error", fontsize=11, weight="bold")
    for row_index, depth in enumerate((1, 4, 10, 20, 30)):
        y_coord = 0.66 - (row_index * 0.14)
        radical = radicals[depth - 1]
        g_partial = g_partials[depth - 1]
        error = abs(g_partial - CATALAN)
        ax_table.text(0.02, y_coord, f"{depth:>2}", family="monospace")
        ax_table.text(0.16, y_coord, f"{radical:.8f}", family="monospace")
        ax_table.text(0.46, y_coord, f"{g_partial:.8f}", family="monospace")
        ax_table.text(0.78, y_coord, f"{error:.2e}", family="monospace")

    _save_figure(fig, path, dpi=180)


CATALAN_CONSTANTS = [
    ("G", CATALAN),
    ("log 2", math.log(2.0)),
    ("G log 2", CATALAN * math.log(2.0)),
    ("G / 2", CATALAN / 2.0),
    ("(2G/pi) log 2", (2.0 * CATALAN / math.pi) * math.log(2.0)),
]


def _lambda_curve(multiplier: np.ndarray, ceiling: float, rate: float) -> np.ndarray:
    return ceiling * (1.0 - np.exp(-rate * (multiplier - 1.0)))


def _plot_flight_map(path: Path, title: str, ceiling: float, rate: float, transposed: bool) -> None:
    multiplier = np.linspace(1.0, 7.0, 400)
    spectrum = _lambda_curve(multiplier, ceiling, rate)

    fig, ax = plt.subplots(figsize=(7.4, 4.7), constrained_layout=True)
    if transposed:
        ax.plot(spectrum, multiplier, color="#155a9c", linewidth=2.2)
        for label, constant in CATALAN_CONSTANTS:
            color = "#c64040" if constant > ceiling else "#239b56"
            ax.axvline(constant, color=color, linestyle="--", linewidth=1.2, alpha=0.85)
            ax.text(constant + 0.01, 6.75, label, rotation=90, va="top", fontsize=8, color=color)
        ax.set_xlabel("Lyapunov lambda")
        ax.set_ylabel("multiplier m1")
    else:
        ax.plot(multiplier, spectrum, color="#155a9c", linewidth=2.2)
        for label, constant in CATALAN_CONSTANTS:
            color = "#c64040" if constant > ceiling else "#239b56"
            ax.axhline(constant, color=color, linestyle="--", linewidth=1.2, alpha=0.85)
            ax.text(6.95, constant + 0.008, label, fontsize=8, color=color, ha="right")
        ax.set_xlabel("multiplier m1")
        ax.set_ylabel("Lyapunov lambda")

    ax.set_title(title)
    ax.grid(alpha=0.22)
    ax.set_xlim(left=0.0 if transposed else 1.0)
    ax.set_ylim(bottom=0.0)

    _save_figure(fig, path, dpi=180)


def _generate_mg_artifacts(paths, records: list[dict[str, object]]) -> None:
    discrete_rows = build_discrete_lattice(2, 12)
    continuous_rows = build_continuous_lattice(2, 12, samples=320)

    discrete_csv = paths.data / "mg_lattice_discrete.csv"
    continuous_csv = paths.data / "mg_lattice_continuous.csv"
    export_csv(discrete_rows, discrete_csv)
    export_csv(continuous_rows, continuous_csv)
    _record(records, "data", discrete_csv.relative_to(paths.root))
    _record(records, "data", continuous_csv.relative_to(paths.root))

    file_plots = [
        ("mg_lattice_checkerboard.png", lambda target: _plot_checkerboard(discrete_rows, target)),
        ("mg_lattice_double_helix.png", lambda target: _plot_double_helix(continuous_rows, target)),
        (
            "mg_lattice_prime_style_funnel.png",
            lambda target: _plot_prime_style_funnel(
                continuous_rows, target, r"$\mathcal{M}_G$ Prime-Style Twisted Funnel Projection"
            ),
        ),
        ("catalan_apery_3d_heatmap.png", lambda target: _plot_catalan_apery_heatmap(discrete_rows, target)),
        ("catalan_Mg_3d.png", lambda target: _plot_dark_mg_matrix(discrete_rows, target)),
        (
            "Mg_catalan_zeta.png",
            lambda target: _plot_prime_style_funnel(
                continuous_rows, target, r"$\mathcal{M}_G$ Prime-Style Twisted Funnel Projection"
            ),
        ),
    ]

    for filename, renderer in file_plots:
        target = paths.figures / filename
        renderer(target)
        _record(records, "figure", target.relative_to(paths.root))


def _generate_prime_artifacts(paths, records: list[dict[str, object]]) -> None:
    figure_1 = paths.figures / "Figure_1.png"
    tanh_path = paths.figures / "tanh_diophantine.png"
    ti2_path = paths.figures / "ti2_primes.png"

    _plot_prime_projection(figure_1)
    _plot_tanh_diophantine(tanh_path)
    _plot_ti2_primes(ti2_path)

    for target in (figure_1, tanh_path, ti2_path):
        _record(records, "figure", target.relative_to(paths.root))


def _generate_pressure_artifacts(paths, records: list[dict[str, object]]) -> None:
    phase_path = paths.figures / "pressure_phase_diagram.png"
    slice_path = paths.figures / "pressure_slice_t1.png"
    _plot_pressure_phase_diagram(phase_path)
    _plot_pressure_slice(slice_path)
    _record(records, "figure", phase_path.relative_to(paths.root))
    _record(records, "figure", slice_path.relative_to(paths.root))


def _generate_catalan_radical_artifacts(paths, records: list[dict[str, object]]) -> None:
    board_png = paths.figures / "catalan_radical_riemann.png"
    spiral_svg = paths.figures / "ramanujan_catalan_spiral_projection.svg"
    _plot_catalan_radical_riemann(board_png)
    _plot_ramanujan_spiral_svg(spiral_svg)
    _record(records, "figure", board_png.relative_to(paths.root))
    _record(records, "figure", spiral_svg.relative_to(paths.root))


def _generate_zded_artifacts(paths, records: list[dict[str, object]]) -> None:
    figure_path = paths.figures / "zded_mg_lattice_calibration.png"
    csv_path = paths.data / "zded_irrational_hypotheses.csv"
    _plot_zded_calibration(figure_path, csv_path)
    _record(records, "figure", figure_path.relative_to(paths.root))
    _record(records, "data", csv_path.relative_to(paths.root))


def _generate_html_svgs(paths, records: list[dict[str, object]]) -> None:
    variants = [
        ("flight_map_lyapunov_vs_catalan.svg", "Lyapunov spectrum vs Catalan constants", math.log(2.0), 0.85, False),
        ("flight_map_lyapunov_vs_catalan_bijective.svg", "Bijective coding (full Bernoulli shift)", math.log(2.0), 0.95, False),
        ("flight_map_lyapunov_vs_catalan_bijective_transposed.svg", "Bijective coding (transposed view)", math.log(2.0), 0.95, True),
        ("flight_map_lyapunov_vs_catalan_injective.svg", "Injective coding (golden-mean ceiling)", math.log(GOLDEN_RATIO), 1.05, False),
        ("flight_map_lyapunov_vs_catalan_injective_transposed.svg", "Injective coding (transposed view)", math.log(GOLDEN_RATIO), 1.05, True),
        ("flight_map_lyapunov_vs_catalan_surjective.svg", "Surjective overlap regime", math.log(3.0), 0.7, False),
    ]

    for filename, title, ceiling, rate, transposed in variants:
        target = paths.html / filename
        _plot_flight_map(target, title, ceiling, rate, transposed)
        _record(records, "html", target.relative_to(paths.root))


SECTION_HANDLERS = {
    "mg": _generate_mg_artifacts,
    "prime": _generate_prime_artifacts,
    "pressure": _generate_pressure_artifacts,
    "catalan": _generate_catalan_radical_artifacts,
    "zded": _generate_zded_artifacts,
    "html": _generate_html_svgs,
}


def generate_all_artifacts(
    output_root: Path | str = "results",
    layout: str = "canonical",
    sections: list[str] | None = None,
) -> list[dict[str, object]]:
    root = Path(output_root)
    paths = resolve_output_paths(root, layout)
    requested_sections = sections or list(SECTION_HANDLERS.keys())

    records: list[dict[str, object]] = []
    for section in requested_sections:
        SECTION_HANDLERS[section](paths, records)

    manifest_path = paths.root / "generated_manifest.csv"
    write_manifest(records, manifest_path)
    _record(records, "meta", manifest_path.relative_to(paths.root))
    write_manifest(records, manifest_path)
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate Catalan Mg figures and data artifacts.")
    parser.add_argument(
        "--output-root",
        default="results",
        help="Target directory. Use '.' with --layout legacy to mirror the historical layout.",
    )
    parser.add_argument(
        "--layout",
        choices=("canonical", "legacy"),
        default="canonical",
        help="canonical writes figures/data/html into subdirectories; legacy writes into the repo-root style layout.",
    )
    parser.add_argument(
        "--section",
        action="append",
        choices=sorted(SECTION_HANDLERS.keys()),
        help="Generate only the named section. Repeat for multiple sections.",
    )

    args = parser.parse_args(argv)
    generate_all_artifacts(args.output_root, args.layout, args.section)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
