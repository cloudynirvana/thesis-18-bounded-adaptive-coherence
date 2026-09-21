#!/usr/bin/env python3
"""Toy coupling tensor and a grounded λ_min margin.

The trajectories are closed-form weight paths on a five-index matrix.
They are not measurements of aging or cancer. Seed 20260921 is used only
for the monotonicity draws.

Research only. Not a medical device.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921

LABELS = [
    "molecular",
    "cellular",
    "tissue",
    "organism",
    "evolutionary",
]
N = 5
ORG = 3
GROUND = 4  # evolutionary index; fixed before any classification

# Predeclared schedules. Not fitted.
DELTA_AGE = 0.04
DELTA_CUT = 0.08
DELTA_AMB_CUT = 0.06
DELTA_AMB_BLOCK = 0.02
DELTA_DIAG = 0.05
SIGMA_HOLD = 0.15
SIGMA_RISE_INTERCEPT = 0.05
SIGMA_RISE_SLOPE = 0.003
T_END = 80.0
N_TIME = 801
DIAG0 = 0.90
CELL_DIAG_PAINTED = 0.95

# Classifier thresholds, fixed before the trajectories are labelled.
THR_LAM_FALL = 0.5
THR_BOTH_FALL = 0.5
THR_SECTOR_BAND = 1.5
THR_CUT_COLLAPSE = 0.25
THR_BLOCK_HELD = 0.9
THR_SECTOR_RISE = 3.0


def initial_weights() -> np.ndarray:
    """Symmetric off-diagonal weights. Diagonal is unused by the Laplacian."""
    pairs = {
        (0, 1): 0.80,
        (0, 2): 0.35,
        (0, 3): 0.55,
        (0, 4): 0.20,
        (1, 2): 0.75,
        (1, 3): 0.60,
        (1, 4): 0.25,
        (2, 3): 0.70,
        (2, 4): 0.30,
        (3, 4): 0.40,
    }
    w = np.zeros((N, N), dtype=float)
    for (i, j), value in pairs.items():
        w[i, j] = w[j, i] = value
    return w


def asymmetric_coupling(weights: np.ndarray) -> np.ndarray:
    """A directional perturbation of W. The Laplacian does not see it."""
    skew = [(1, 3, 0.04), (2, 3, -0.03), (0, 1, 0.02)]
    c = weights.copy()
    for i, j, amp in skew:
        c[i, j] = weights[i, j] + amp
        c[j, i] = weights[j, i] - amp
    np.fill_diagonal(c, DIAG0)
    return c


def symmetrize(weights: np.ndarray) -> np.ndarray:
    a = np.array(weights, dtype=float, copy=True)
    a = 0.5 * (a + a.T)
    np.fill_diagonal(a, 0.0)
    return a


def laplacian(weights: np.ndarray) -> np.ndarray:
    a = symmetrize(weights)
    degree = a.sum(axis=1)
    return np.diag(degree) - a


def algebraic_connectivity(weights: np.ndarray) -> float:
    ev = np.linalg.eigvalsh(laplacian(weights))
    return float(np.sort(ev)[1])


def fiedler_stats(weights: np.ndarray) -> tuple[float, float]:
    ev, vec = np.linalg.eigh(laplacian(weights))
    order = np.argsort(ev)
    v = vec[:, order[1]]
    v = v / np.linalg.norm(v)
    participation = float(1.0 / np.sum(v**4))
    organism_mass = float(v[ORG] ** 2)
    return participation, organism_mass


def grounded_lambda_min(weights: np.ndarray, ground: int = GROUND) -> float:
    ell = laplacian(weights)
    keep = [i for i in range(N) if i != ground]
    gamma = ell[np.ix_(keep, keep)]
    return float(np.linalg.eigvalsh(gamma)[0])


def induced_block_lambda2(weights: np.ndarray, drop: int = ORG) -> float:
    """Algebraic connectivity of the induced subgraph on all indices but `drop`."""
    keep = [i for i in range(N) if i != drop]
    block = symmetrize(weights)[np.ix_(keep, keep)]
    ev = np.linalg.eigvalsh(laplacian(block))
    return float(np.sort(ev)[1])


def cut_and_block_means(weights: np.ndarray) -> tuple[float, float, float]:
    a = symmetrize(weights)
    cut = []
    block = []
    for i in range(N):
        for j in range(i + 1, N):
            if i == ORG or j == ORG:
                cut.append(a[i, j])
            else:
                block.append(a[i, j])
    cut_sum = float(np.sum(cut))
    return float(np.mean(cut)), float(np.mean(block)), cut_sum


def edge_rates(kind: str) -> np.ndarray:
    rates = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i + 1, N):
            incident = i == ORG or j == ORG
            if kind == "held":
                rate = 0.0
            elif kind == "aging":
                rate = DELTA_AGE
            elif kind == "cancer":
                rate = DELTA_CUT if incident else 0.0
            elif kind == "ambiguous":
                rate = DELTA_AMB_CUT if incident else DELTA_AMB_BLOCK
            elif kind == "diagonal":
                rate = 0.0
            else:
                raise ValueError(kind)
            rates[i, j] = rates[j, i] = rate
    return rates


def weights_at(t: float, rates: np.ndarray, w0: np.ndarray) -> np.ndarray:
    w = np.zeros_like(w0)
    for i in range(N):
        for j in range(i + 1, N):
            value = w0[i, j] * np.exp(-rates[i, j] * t)
            w[i, j] = w[j, i] = value
    return w


def diagonal_at(t: float, kind: str) -> np.ndarray:
    d = np.full(N, DIAG0, dtype=float)
    if kind == "aging":
        d *= np.exp(-DELTA_DIAG * t)
    elif kind == "diagonal":
        d *= np.exp(-DELTA_DIAG * t)
    elif kind == "cancer":
        d[1] = CELL_DIAG_PAINTED
    return d


def rayleigh(weights: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    a = symmetrize(weights)
    quad = float(x @ laplacian(weights) @ x)
    acc = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            acc += a[i, j] * (x[i] - x[j]) ** 2
    return quad, float(acc)


def classify(lam_ratio: float, block_ratio: float, cut_ratio: float, sector_change: float) -> str:
    aging = (
        lam_ratio < THR_LAM_FALL
        and block_ratio < THR_BOTH_FALL
        and cut_ratio < THR_BOTH_FALL
        and (1.0 / THR_SECTOR_BAND) < sector_change < THR_SECTOR_BAND
    )
    cancer = (
        lam_ratio < THR_LAM_FALL
        and cut_ratio < THR_CUT_COLLAPSE
        and block_ratio > THR_BLOCK_HELD
        and sector_change > THR_SECTOR_RISE
    )
    if aging and cancer:
        return "both"
    if aging:
        return "aging-like"
    if cancer:
        return "cancer-like"
    return "neither"


def first_crossing(t: np.ndarray, lam: np.ndarray, sigma: np.ndarray) -> float | None:
    v = lam - sigma
    if np.all(v > 0):
        return None
    idx = int(np.argmax(v <= 0))
    if idx == 0:
        return float(t[0])
    t0, t1 = float(t[idx - 1]), float(t[idx])
    v0, v1 = float(v[idx - 1]), float(v[idx])
    if v1 == v0:
        return t1
    return t0 + (0.0 - v0) * (t1 - t0) / (v1 - v0)


def rho_family(n_rho: int = 20) -> dict:
    rhos = np.linspace(0.0, 0.95, n_rho)
    raw = []
    alg = []
    grounded = []
    for rho in rhos:
        c = (1.0 - rho) * np.eye(N) + rho * np.ones((N, N))
        raw.append(float(np.sort(np.linalg.eigvalsh(0.5 * (c + c.T)))[0]))
        # Off-diagonal coupling in this family is ρ. Diagonal is irrelevant to L.
        alg.append(algebraic_connectivity(c))
        grounded.append(grounded_lambda_min(c))
    raw_a = np.array(raw)
    alg_a = np.array(alg)
    return {
        "rho": [float(x) for x in rhos],
        "lambda_min_raw": raw,
        "lambda2": alg,
        "lambda_min_grounded": grounded,
        "raw_decreases": bool(np.all(np.diff(raw_a) < 0)),
        "lambda2_increases": bool(np.all(np.diff(alg_a) > 0)),
        "analytic_raw_max_abs_err": float(np.max(np.abs(raw_a - (1.0 - rhos)))),
        "analytic_lambda2_max_abs_err": float(np.max(np.abs(alg_a - (N * rhos)))),
    }


def cycle_counterexample() -> dict:
    r = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    ev = np.linalg.eigvals(r)
    return {
        "matrix": r.tolist(),
        "eigenvalues_real": [float(z.real) for z in ev],
        "eigenvalues_imag": [float(z.imag) for z in ev],
        "has_nonreal": bool(np.any(np.abs(ev.imag) > 1e-10)),
    }


def run_kind(kind: str, t: np.ndarray, w0: np.ndarray) -> dict:
    rates = edge_rates(kind)
    lam2 = np.empty_like(t)
    lmin = np.empty_like(t)
    lblock = np.empty_like(t)
    mean_cut = np.empty_like(t)
    mean_block = np.empty_like(t)
    cut_sum = np.empty_like(t)
    part = np.empty_like(t)
    mass = np.empty_like(t)
    dmin = np.empty_like(t)
    for k, time in enumerate(t):
        w = weights_at(float(time), rates, w0)
        lam2[k] = algebraic_connectivity(w)
        lmin[k] = grounded_lambda_min(w, GROUND)
        lblock[k] = induced_block_lambda2(w, ORG)
        mean_cut[k], mean_block[k], cut_sum[k] = cut_and_block_means(w)
        part[k], mass[k] = fiedler_stats(w)
        dmin[k] = float(np.min(diagonal_at(float(time), kind)))
    sector0 = mean_block[0] / mean_cut[0]
    sector_t = mean_block[-1] / mean_cut[-1]
    label = classify(
        float(lmin[-1] / lmin[0]),
        float(mean_block[-1] / mean_block[0]),
        float(mean_cut[-1] / mean_cut[0]),
        float(sector_t / sector0),
    )
    sigma_hold = np.full_like(t, SIGMA_HOLD)
    sigma_rise = SIGMA_RISE_INTERCEPT + SIGMA_RISE_SLOPE * t
    # Cut bound for λ2: n / (n_A n_B) * Cut with n_A = 1.
    bound = (N / (1 * (N - 1))) * cut_sum
    # Grounded test: λ_min(Γ) ≤ Cut, ground ≠ organism.
    return {
        "kind": kind,
        "class": label,
        "t": t,
        "lambda2": lam2,
        "lambda_min": lmin,
        "lambda2_block": lblock,
        "mean_cut": mean_cut,
        "mean_block": mean_block,
        "cut_sum": cut_sum,
        "participation": part,
        "organism_mass": mass,
        "min_diagonal": dmin,
        "sigma_hold": sigma_hold,
        "sigma_rise": sigma_rise,
        "crossing_hold": first_crossing(t, lmin, sigma_hold),
        "crossing_rise": first_crossing(t, lmin, sigma_rise),
        "lambda2_bound_slack_min": float(np.min(bound - lam2)),
        "grounded_cut_slack_min": float(np.min(cut_sum - lmin)),
        "lambda_min_ratio": float(lmin[-1] / lmin[0]),
        "lambda2_block_ratio": float(lblock[-1] / lblock[0]),
        "mean_cut_ratio": float(mean_cut[-1] / mean_cut[0]),
        "mean_block_ratio": float(mean_block[-1] / mean_block[0]),
        "sector_ratio_0": float(sector0),
        "sector_ratio_T": float(sector_t),
        "sector_change": float(sector_t / sector0),
    }


def snapshot_rows(series: dict, times: list[float]) -> list[dict]:
    rows = []
    t = series["t"]
    for time in times:
        k = int(np.argmin(np.abs(t - time)))
        rows.append(
            {
                "t": float(t[k]),
                "lambda_min": float(series["lambda_min"][k]),
                "lambda2": float(series["lambda2"][k]),
                "lambda2_block": float(series["lambda2_block"][k]),
                "mean_cut": float(series["mean_cut"][k]),
                "mean_block": float(series["mean_block"][k]),
                "organism_mass": float(series["organism_mass"][k]),
                "participation": float(series["participation"][k]),
                "min_diagonal": float(series["min_diagonal"][k]),
                "V_hold": float(series["lambda_min"][k] - SIGMA_HOLD),
            }
        )
    return rows


def checks(w0: np.ndarray, runs: dict[str, dict], family: dict, cycle: dict) -> dict:
    rng = np.random.default_rng(SEED)
    mono_fail = 0
    grounded_fail = 0
    for _ in range(200):
        i, j = rng.choice(N, size=2, replace=False)
        lifted = w0.copy()
        lifted[i, j] += 0.05
        lifted[j, i] += 0.05
        if algebraic_connectivity(lifted) + 1e-10 < algebraic_connectivity(w0):
            mono_fail += 1
        if grounded_lambda_min(lifted) + 1e-10 < grounded_lambda_min(w0):
            grounded_fail += 1

    ray_err = 0.0
    for _ in range(50):
        x = rng.normal(size=N)
        x = x - x.mean()
        x = x / np.linalg.norm(x)
        quad, acc = rayleigh(w0, x)
        ray_err = max(ray_err, abs(quad - acc))

    aging = runs["aging"]
    scale = aging["lambda_min"][0] * np.exp(-DELTA_AGE * aging["t"])
    scale2 = aging["lambda2"][0] * np.exp(-DELTA_AGE * aging["t"])
    scale_b = aging["lambda2_block"][0] * np.exp(-DELTA_AGE * aging["t"])
    cancer = runs["cancer"]
    # Painting the cellular diagonal must not move the grounded eigenvalue.
    w_paint = w0.copy()
    # diagonals are already zero in W; compare two C diagonals through the same W.
    paint_gap = abs(grounded_lambda_min(w_paint) - cancer["lambda_min"][0])

    c_asym = asymmetric_coupling(w0)
    ev_asym = np.linalg.eigvals(c_asym)
    # Laplacian of C equals Laplacian of its symmetric off-diagonal part.
    gap_asym = abs(grounded_lambda_min(c_asym) - grounded_lambda_min(w0))

    expected = {
        "held": "neither",
        "aging": "aging-like",
        "cancer": "cancer-like",
        "ambiguous": "neither",
        "diagonal": "neither",
    }
    class_ok = {name: runs[name]["class"] == label for name, label in expected.items()}

    analytic_cross = float(np.log(aging["lambda_min"][0] / SIGMA_HOLD) / DELTA_AGE)
    cross_err = abs(aging["crossing_hold"] - analytic_cross)

    return {
        "monotonicity_lambda2_failures": mono_fail,
        "monotonicity_grounded_failures": grounded_fail,
        "rayleigh_max_abs_err": float(ray_err),
        "aging_lambda_min_scale_max_abs_err": float(np.max(np.abs(aging["lambda_min"] - scale))),
        "aging_lambda2_scale_max_abs_err": float(np.max(np.abs(aging["lambda2"] - scale2))),
        "aging_block_scale_max_abs_err": float(np.max(np.abs(aging["lambda2_block"] - scale_b))),
        "diagonal_lambda_min_range": float(np.max(runs["diagonal"]["lambda_min"]) - np.min(runs["diagonal"]["lambda_min"])),
        "cancer_block_lambda2_range": float(np.max(cancer["lambda2_block"]) - np.min(cancer["lambda2_block"])),
        "painted_diagonal_gap_at_0": float(paint_gap),
        "asymmetry_grounded_gap": float(gap_asym),
        "asymmetry_has_nonreal": bool(np.any(np.abs(ev_asym.imag) > 1e-8)),
        "asymmetry_eigenvalues_real": [float(z.real) for z in ev_asym],
        "asymmetry_eigenvalues_imag": [float(z.imag) for z in ev_asym],
        "cycle_has_nonreal": cycle["has_nonreal"],
        "rho_raw_decreases": family["raw_decreases"],
        "rho_lambda2_increases": family["lambda2_increases"],
        "rho_raw_analytic_err": family["analytic_raw_max_abs_err"],
        "rho_lambda2_analytic_err": family["analytic_lambda2_max_abs_err"],
        "classes_match_predeclaration": class_ok,
        "aging_analytic_crossing": analytic_cross,
        "aging_numeric_crossing_abs_err": float(cross_err),
        "cancer_lambda2_bound_slack_min": cancer["lambda2_bound_slack_min"],
        "cancer_grounded_cut_slack_min": cancer["grounded_cut_slack_min"],
        "aging_lambda2_bound_slack_min": aging["lambda2_bound_slack_min"],
        "all_slacks_nonnegative": bool(
            cancer["lambda2_bound_slack_min"] >= -1e-9
            and cancer["grounded_cut_slack_min"] >= -1e-9
            and aging["lambda2_bound_slack_min"] >= -1e-9
        ),
    }


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )


COLORS = {
    "held": "#4C5760",
    "aging": "#1F4E79",
    "cancer": "#8C3A3A",
    "ambiguous": "#8A6A12",
    "diagonal": "#2F6B4F",
}
NAMES = {
    "held": "held weights",
    "aging": "global decay",
    "cancer": "organism-scale cut",
    "ambiguous": "uneven decay",
    "diagonal": "diagonal only",
}


def draw(runs: dict[str, dict], family: dict) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    style()

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for key in ("held", "aging", "cancer", "ambiguous", "diagonal"):
        s = runs[key]
        ax.plot(s["t"], s["lambda_min"], color=COLORS[key], lw=1.6, label=NAMES[key])
    ax.axhline(SIGMA_HOLD, color="#666666", ls="--", lw=0.9, label="constant proxy σ = 0.15")
    ax.set_xlabel("toy time")
    ax.set_ylabel("λ_min of the grounded operator")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "lambda_min_trajectories.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for key, ls_cut, ls_block in (("aging", "-", "--"), ("cancer", "-", "--")):
        s = runs[key]
        ax.plot(s["t"], s["mean_cut"], color=COLORS[key], ls=ls_cut, lw=1.6, label=f"{NAMES[key]}, cut mean")
        ax.plot(s["t"], s["mean_block"], color=COLORS[key], ls=ls_block, lw=1.6, label=f"{NAMES[key]}, block mean")
    ax.set_xlabel("toy time")
    ax.set_ylabel("mean symmetric weight")
    ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "sector_means.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    rho = np.array(family["rho"])
    ax.plot(rho, family["lambda_min_raw"], color="#8C3A3A", lw=1.7, label="λ_min of C(ρ)")
    ax.plot(rho, family["lambda2"], color="#1F4E79", lw=1.7, label="λ₂ of the Laplacian")
    ax.set_xlabel("coupling parameter ρ")
    ax.set_ylabel("eigenvalue")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "raw_eigenvalue_versus_coupling.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for key in ("aging", "cancer"):
        s = runs[key]
        ax.plot(s["t"], s["organism_mass"], color=COLORS[key], lw=1.6, label=f"{NAMES[key]}, organism mass")
    ax.set_xlabel("toy time")
    ax.set_ylabel("squared Fiedler component on the organism index")
    ax.set_ylim(0.0, 1.0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fiedler_organism_mass.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    aging = runs["aging"]
    cancer = runs["cancer"]
    ax.plot(aging["t"], aging["lambda_min"] - SIGMA_HOLD, color=COLORS["aging"], lw=1.6, label="global decay, σ = 0.15")
    ax.plot(
        aging["t"],
        aging["lambda_min"] - aging["sigma_rise"],
        color=COLORS["aging"],
        ls="--",
        lw=1.4,
        label="global decay, rising proxy",
    )
    ax.plot(cancer["t"], cancer["lambda_min"] - SIGMA_HOLD, color=COLORS["cancer"], lw=1.6, label="organism-scale cut, σ = 0.15")
    ax.axhline(0.0, color="#444444", lw=0.8)
    ax.set_xlabel("toy time")
    ax.set_ylabel("margin V")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "margin_versus_proxy.png")
    plt.close(fig)


def round_floats(obj, nd=6):
    if isinstance(obj, float):
        return None if obj != obj else round(obj, nd)
    if isinstance(obj, dict):
        return {k: round_floats(v, nd) for k, v in obj.items()}
    if isinstance(obj, list):
        return [round_floats(v, nd) for v in obj]
    if isinstance(obj, (np.floating,)):
        return round(float(obj), nd)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def main() -> None:
    w0 = initial_weights()
    t = np.linspace(0.0, T_END, N_TIME)
    kinds = ("held", "aging", "cancer", "ambiguous", "diagonal")
    runs = {kind: run_kind(kind, t, w0) for kind in kinds}
    family = rho_family()
    cycle = cycle_counterexample()
    report = checks(w0, runs, family, cycle)

    hard = [
        report["monotonicity_lambda2_failures"] == 0,
        report["monotonicity_grounded_failures"] == 0,
        report["rayleigh_max_abs_err"] < 1e-9,
        report["aging_lambda_min_scale_max_abs_err"] < 1e-8,
        report["aging_lambda2_scale_max_abs_err"] < 1e-8,
        report["aging_block_scale_max_abs_err"] < 1e-8,
        report["diagonal_lambda_min_range"] < 1e-10,
        report["cancer_block_lambda2_range"] < 1e-10,
        report["painted_diagonal_gap_at_0"] < 1e-12,
        report["asymmetry_grounded_gap"] < 1e-12,
        report["cycle_has_nonreal"] is True,
        report["rho_raw_decreases"] is True,
        report["rho_lambda2_increases"] is True,
        report["rho_raw_analytic_err"] < 1e-8,
        report["rho_lambda2_analytic_err"] < 1e-8,
        all(report["classes_match_predeclaration"].values()),
        report["aging_numeric_crossing_abs_err"] < 1e-3,
        report["all_slacks_nonnegative"] is True,
        runs["held"]["crossing_hold"] is None,
        runs["diagonal"]["crossing_hold"] is None,
        runs["aging"]["crossing_hold"] is not None,
        runs["cancer"]["crossing_hold"] is not None,
    ]
    if not all(hard):
        raise SystemExit(f"toy checks failed: {json.dumps(round_floats(report), indent=2)}")

    times = [0.0, 20.0, 40.0, 60.0, 80.0]
    c_asym = asymmetric_coupling(w0)
    payload = {
        "seed": SEED,
        "labels": LABELS,
        "organism_index": ORG,
        "ground_index": GROUND,
        "ground_label": LABELS[GROUND],
        "W0": w0.tolist(),
        "C0_asymmetric": c_asym.tolist(),
        "parameters": {
            "delta_age": DELTA_AGE,
            "delta_cut": DELTA_CUT,
            "delta_ambiguous_cut": DELTA_AMB_CUT,
            "delta_ambiguous_block": DELTA_AMB_BLOCK,
            "delta_diagonal": DELTA_DIAG,
            "sigma_hold": SIGMA_HOLD,
            "sigma_rise_intercept": SIGMA_RISE_INTERCEPT,
            "sigma_rise_slope": SIGMA_RISE_SLOPE,
            "t_end": T_END,
            "n_time": N_TIME,
            "diag0": DIAG0,
            "cellular_diagonal_painted": CELL_DIAG_PAINTED,
        },
        "thresholds": {
            "lambda_ratio_fall": THR_LAM_FALL,
            "both_sectors_fall": THR_BOTH_FALL,
            "sector_change_band": THR_SECTOR_BAND,
            "cut_collapse": THR_CUT_COLLAPSE,
            "block_held": THR_BLOCK_HELD,
            "sector_rise": THR_SECTOR_RISE,
        },
        "initial": {
            "lambda_min": float(runs["held"]["lambda_min"][0]),
            "lambda2": float(runs["held"]["lambda2"][0]),
            "lambda2_block": float(runs["held"]["lambda2_block"][0]),
            "mean_cut": float(runs["held"]["mean_cut"][0]),
            "mean_block": float(runs["held"]["mean_block"][0]),
            "cut_sum": float(runs["held"]["cut_sum"][0]),
            "organism_mass": float(runs["held"]["organism_mass"][0]),
            "participation": float(runs["held"]["participation"][0]),
        },
        "summaries": {},
        "snapshots": {},
        "rho_endpoints": {
            "rho0_raw": family["lambda_min_raw"][0],
            "rho0_lambda2": family["lambda2"][0],
            "rho095_raw": family["lambda_min_raw"][-1],
            "rho095_lambda2": family["lambda2"][-1],
        },
        "cycle": cycle,
        "checks": report,
    }
    for kind, series in runs.items():
        payload["summaries"][kind] = {
            "class": series["class"],
            "crossing_hold": series["crossing_hold"],
            "crossing_rise": series["crossing_rise"],
            "lambda_min_ratio": series["lambda_min_ratio"],
            "lambda2_block_ratio": series["lambda2_block_ratio"],
            "mean_cut_ratio": series["mean_cut_ratio"],
            "mean_block_ratio": series["mean_block_ratio"],
            "sector_ratio_0": series["sector_ratio_0"],
            "sector_ratio_T": series["sector_ratio_T"],
            "sector_change": series["sector_change"],
            "lambda2_bound_slack_min": series["lambda2_bound_slack_min"],
            "grounded_cut_slack_min": series["grounded_cut_slack_min"],
            "organism_mass_0": float(series["organism_mass"][0]),
            "organism_mass_T": float(series["organism_mass"][-1]),
            "participation_0": float(series["participation"][0]),
            "participation_T": float(series["participation"][-1]),
            "min_diagonal_T": float(series["min_diagonal"][-1]),
            "lambda_min_T": float(series["lambda_min"][-1]),
            "lambda2_T": float(series["lambda2"][-1]),
            "lambda2_block_T": float(series["lambda2_block"][-1]),
        }
        payload["snapshots"][kind] = snapshot_rows(series, times)

    draw(runs, family)
    out = ROOT / "results.json"
    out.write_text(json.dumps(round_floats(payload), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    for kind in kinds:
        s = payload["summaries"][kind]
        print(
            f"{kind:10} {s['class']:12} cross={s['crossing_hold']} "
            f"lmin_ratio={s['lambda_min_ratio']:.6f} block_l2_ratio={s['lambda2_block_ratio']:.6f} "
            f"sector_change={s['sector_change']:.4f}"
        )


if __name__ == "__main__":
    main()
