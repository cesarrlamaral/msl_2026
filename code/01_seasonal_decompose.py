"""
Classical time-series decomposition of the MSLRAD TOTALSERIES dose-rate record
(detectors B and E, SOL 00000-02358), replacing the R decompose() pass that
used frequency=100 observations (~1 sol) -- far too short to isolate a real
Martian seasonal (annual, ~668.6 sol) cycle. Here the period is set from the
actual sampling rate of the data (obs/sol), giving a period in observation
units that corresponds to one Martian year.

Implemented by hand (centered moving-average trend + averaged-by-phase
seasonal component + residual) rather than via statsmodels' STL/
seasonal_decompose: both use LOESS/direct convolution windows that scale
with the period, and at period ~67600 (vs. ~238500 total observations)
this was computationally infeasible in practice (two attempts killed after
several minutes with no result). The classical additive decomposition
below is exactly the algorithm R's own decompose() uses, and is O(n) here
via numpy cumulative-sum moving averages, so it completes in seconds.

Author: Cesar R. L. Amaral (analysis re-implemented in Python, 2026-09-05,
per project decision to move off R -- see feedback_mslrad_python_nao_r).
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\MSLRADv1 TOTALSERIES RESULTS"
OUT_DIR = r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\python_reanalysis"

SOL_SPAN = 2358  # SOL 00000-02358, per README TOTALSERIES.txt
MARS_YEAR_SOLS = 668.6  # length of one Martian year in sols


def load_series(detector):
    path = f"{DATA_DIR}\\MSLRAD_TOTALSERIES_{detector}.txt"
    s = pd.read_csv(path, header=0).iloc[:, 0].values.astype(float)
    return s


def centered_moving_average(y, window):
    """Centered moving average, matching R decompose()'s trend filter
    (equal-weight for odd window, half-weight at the two end points for
    even window). Implemented via cumulative sum for O(n) performance."""
    n = len(y)
    half = window // 2
    csum = np.cumsum(np.insert(y, 0, 0))
    trend = np.full(n, np.nan)
    if window % 2 == 1:
        # odd window: simple centered average
        trend[half:n - half] = (csum[window:] - csum[:-window]) / window
    else:
        # even window: R's decompose() uses a 2xwindow filter (half weight at ends)
        w2 = window + 1
        avg = (csum[window:] - csum[:-window]) / window
        trend[half:n - half] = (avg[:-1] + avg[1:]) / 2
    return trend


def classical_decompose(y, period):
    n = len(y)
    trend = centered_moving_average(y, period)
    detrended = y - trend

    # Average detrended value at each phase (index mod period) across all
    # available cycles, ignoring NaNs at the series edges (no trend there).
    phase = np.arange(n) % period
    seasonal_avg = np.full(period, np.nan)
    for p in range(period):
        vals = detrended[phase == p]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            seasonal_avg[p] = vals.mean()
    # center the seasonal component around zero, as decompose() does
    seasonal_avg -= np.nanmean(seasonal_avg)
    seasonal = seasonal_avg[phase]

    resid = y - trend - seasonal
    return trend, seasonal, resid


def decompose_and_plot(detector):
    y = load_series(detector)
    n = len(y)
    obs_per_sol = n / SOL_SPAN
    period = round(MARS_YEAR_SOLS * obs_per_sol)
    n_cycles = n / period

    print(f"[{detector}] n={n} obs, obs/sol={obs_per_sol:.3f}, "
          f"annual period={period} obs, cycles available={n_cycles:.2f}")

    trend, seasonal, resid = classical_decompose(y, period)
    sol_axis = np.arange(n) / obs_per_sol

    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(sol_axis, y, lw=0.3, color="navy")
    axes[0].set_ylabel("Observed\n(uGy/hour)")
    axes[0].set_title(f"Classical additive decomposition, detector {detector} "
                       f"(period = {period} obs \u2248 1 Martian year)")

    axes[1].plot(sol_axis, trend, lw=1, color="darkred")
    axes[1].set_ylabel("Trend")

    axes[2].plot(sol_axis, seasonal, lw=0.5, color="darkgreen")
    axes[2].set_ylabel("Seasonal\n(annual)")

    axes[3].plot(sol_axis, resid, lw=0.3, color="gray")
    axes[3].set_ylabel("Residual")
    axes[3].set_xlabel("SOL (mission day)")

    fig.tight_layout()
    out_png = f"{OUT_DIR}\\decomposition_{detector}_annual_period.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    out_df = pd.DataFrame({
        "sol": sol_axis, "observed": y, "trend": trend,
        "seasonal": seasonal, "resid": resid,
    })
    out_csv = f"{OUT_DIR}\\decomposition_{detector}_annual_period.csv"
    out_df.to_csv(out_csv, index=False)
    print(f"  saved {out_csv}")

    seasonal_amp = np.nanstd(seasonal)
    resid_amp = np.nanstd(resid)
    print(f"  seasonal component std = {seasonal_amp:.4f}, residual std = {resid_amp:.4f}, "
          f"ratio = {seasonal_amp/resid_amp:.4f}")

    # Also plot just one cycle of the seasonal component to see its actual shape
    fig2, ax2 = plt.subplots(figsize=(8, 3.5))
    sol_phase = np.arange(period) / obs_per_sol
    seasonal_one_cycle = seasonal[:period]
    ax2.plot(sol_phase, seasonal_one_cycle, color="darkgreen")
    ax2.set_xlabel("Sol within the ~1-Martian-year cycle")
    ax2.set_ylabel("Seasonal component (uGy/hour)")
    ax2.set_title(f"Detector {detector}: shape of the fitted annual-period component")
    fig2.tight_layout()
    out_png2 = f"{OUT_DIR}\\decomposition_{detector}_seasonal_shape.png"
    fig2.savefig(out_png2, dpi=200, bbox_inches="tight")
    plt.close(fig2)
    print(f"  saved {out_png2}")

    return trend, seasonal, resid


if __name__ == "__main__":
    for det in ["B", "E"]:
        decompose_and_plot(det)
