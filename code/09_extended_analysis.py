"""
09_extended_analysis.py

Decomposicao classica + wavelet (Morlet) na serie de dose do RAD agora
estendida a sol 1-4843 (~7.2 anos marcianos, contra 3.5 na versao do paper),
usando sol real como eixo de tempo (nao mais "indice de observacao" como
proxy, ja que agora cada observacao tem sol verdadeiro -- ver
07_backfill_original_history.py).

Substitui/estende 01_seasonal_decompose.py e 02_wavelet_totalseries.py.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pycwt as wavelet

from pathlib import Path

HERE = Path(__file__).parent
MASTER = HERE / "mslrad_master_sol_series.csv"
MARS_YEAR_SOLS = 668.6

# The classical phase-averaged seasonal component is NOT robust to Solar
# Particle Events: a single-sol SPE excursion falls in one phase bin of the
# 669-sol modulo average, and because the reconstructed seasonal series tiles
# that per-phase average across every cycle (seasonal = seasonal_avg[phase]),
# one real transient event at one actual sol reappears as a fake "spike every
# ~669 sols" across the ENTIRE plotted record -- an artifact of the
# reconstruction, not a detection of anything periodic. The original paper
# already flags heavy-SPE windows separately in Table 1 (e.g. "2014 SPE",
# "2017 SPE") rather than folding them into the secular-trend discussion;
# here that same logic is applied automatically (see flag_spe_outlier_sols
# below) before computing the seasonal/trend/residual decomposition. Excluded
# sols are left untouched in the "observed" panel and in the wavelet (which
# correctly localizes each as a short-period transient, not an annual signal).
SPE_MAD_THRESHOLD = 8.0  # sols with |value - rolling_median| > this many MADs are excluded
SPE_ROLLING_WINDOW = 21  # sols, short window so the local trend doesn't itself absorb the spike


def load_master():
    return pd.read_csv(MASTER, parse_dates=["date"])


def interpolate_gaps(y: np.ndarray) -> tuple[np.ndarray, float]:
    """Linear-interpolate internal NaN gaps; report the fraction interpolated.
    Leading/trailing NaNs (before first / after last real obs) are left as NaN
    and trimmed by the caller."""
    s = pd.Series(y)
    first, last = s.first_valid_index(), s.last_valid_index()
    s = s.iloc[first:last + 1]
    frac_missing = s.isna().mean()
    s_interp = s.interpolate(method="linear", limit_direction="both")
    return s_interp.values, frac_missing, first


def flag_spe_outlier_sols(y: np.ndarray, window: int = SPE_ROLLING_WINDOW,
                           mad_threshold: float = SPE_MAD_THRESHOLD) -> np.ndarray:
    """Robust outlier flag: sols departing from a short rolling median by more
    than mad_threshold robust MADs. Returns a boolean mask, True = outlier."""
    s = pd.Series(y)
    roll_med = s.rolling(window, center=True, min_periods=1).median()
    dev = (s - roll_med).abs()
    mad = 1.4826 * dev.median()  # scale factor makes MAD ~ std for normal data
    if mad == 0 or np.isnan(mad):
        return np.zeros(len(y), dtype=bool)
    return (dev > mad_threshold * mad).values


def centered_moving_average(y, window):
    n = len(y)
    half = window // 2
    csum = np.cumsum(np.insert(y, 0, 0))
    trend = np.full(n, np.nan)
    if window % 2 == 1:
        trend[half:n - half] = (csum[window:] - csum[:-window]) / window
    else:
        avg = (csum[window:] - csum[:-window]) / window
        trend[half:n - half] = (avg[:-1] + avg[1:]) / 2
    return trend


def classical_decompose(y, period):
    n = len(y)
    trend = centered_moving_average(y, period)
    detrended = y - trend
    phase = np.arange(n) % period
    seasonal_avg = np.full(period, np.nan)
    for p in range(period):
        vals = detrended[phase == p]
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            seasonal_avg[p] = vals.mean()
    seasonal_avg -= np.nanmean(seasonal_avg)
    seasonal = seasonal_avg[phase]
    resid = y - trend - seasonal
    return trend, seasonal, resid


def run_for_detector(df: pd.DataFrame, detector: str):
    col = f"dose_{detector}_mean"
    y_raw = df[col].values
    y, frac_missing, first_idx = interpolate_gaps(y_raw)
    sol_axis = df["sol"].values[first_idx:first_idx + len(y)]
    n = len(y)
    n_cycles = (sol_axis[-1] - sol_axis[0]) / MARS_YEAR_SOLS
    print(f"\n=== Detector {detector} ===")
    print(f"n={n} sols (sol {sol_axis[0]}..{sol_axis[-1]}), {frac_missing*100:.1f}% interpolated, "
          f"{n_cycles:.2f} Martian years covered (vs 3.53 in the original paper)")

    # --- classical decomposition, period in real sols now ---
    # SPE-excluded input for decomposition only (see note above)
    spe_mask = flag_spe_outlier_sols(y)
    n_spe = spe_mask.sum()
    spe_sols_list = sol_axis[spe_mask]
    print(f"  {n_spe} sol(s) flagged as SPE outliers (>{SPE_MAD_THRESHOLD:.0f} MAD from local median), "
          f"interpolated out before decomposition: {list(spe_sols_list)}")
    y_decomp = np.where(spe_mask, np.nan, y)
    y_decomp = pd.Series(y_decomp).interpolate(limit_direction="both").values

    period = round(MARS_YEAR_SOLS)
    trend, seasonal, resid = classical_decompose(y_decomp, period)
    seasonal_amp = np.nanstd(seasonal)
    resid_amp = np.nanstd(resid)
    print(f"classical decomposition: seasonal std={seasonal_amp:.4f}, residual std={resid_amp:.4f}, "
          f"ratio={seasonal_amp/resid_amp:.4f}")

    fig, axes = plt.subplots(4, 1, figsize=(13, 9), sharex=True)
    axes[0].plot(sol_axis, y, lw=0.4, color="navy")
    axes[0].set_ylabel("Observed\n(uGy/hour)")
    axes[0].set_title(f"Extended classical decomposition, detector {detector}, sol {sol_axis[0]}-{sol_axis[-1]} "
                       f"({n_cycles:.1f} Martian years; period={period} sols)")
    axes[1].plot(sol_axis, trend, lw=1, color="darkred")
    axes[1].set_ylabel("Trend")
    axes[2].plot(sol_axis, seasonal, lw=0.5, color="darkgreen")
    axes[2].set_ylabel("Seasonal\n(annual)")
    axes[3].plot(sol_axis, resid, lw=0.3, color="gray")
    axes[3].set_ylabel("Residual")
    axes[3].set_xlabel("SOL (mission day)")
    fig.tight_layout()
    out_png = HERE / f"decomposition_{detector}_extended_full_history.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    # --- wavelet, dt = 1 sol (real units, no more obs/sol scaling factor) ---
    y0 = y - y.mean()
    y_norm = y0 / y0.std()
    dt = 1.0
    mother = wavelet.Morlet(6)
    s0 = 2 * dt
    dj = 1 / 12
    J = int(np.log2(n / 3 / s0) / dj)
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(y_norm, dt, dj, s0, J, mother)
    power = np.abs(wave) ** 2
    period_sols = 1 / freqs
    alpha, _, _ = wavelet.ar1(y_norm)
    signif, _ = wavelet.significance(1.0, dt, scales, 0, alpha, significance_level=0.95, wavelet=mother)
    sig95 = power / (np.ones([1, n]) * signif[:, None])

    fig2, ax2 = plt.subplots(figsize=(14, 6))
    levels = np.linspace(0, np.percentile(power, 99.5), 30)
    cf = ax2.contourf(sol_axis, np.log2(period_sols), power, levels=levels, extend="both", cmap="viridis")
    ax2.contour(sol_axis, np.log2(period_sols), sig95, [1], colors="white", linewidths=1.2)
    ax2.fill(np.concatenate([sol_axis, sol_axis[-1:], sol_axis[:1]]),
              np.concatenate([np.log2(coi), [np.log2(period_sols.max())], [np.log2(period_sols.max())]]),
              "k", alpha=0.4, hatch="/")
    ax2.axhline(np.log2(MARS_YEAR_SOLS), color="red", ls="--", lw=1,
                label=f"1 Martian year ({MARS_YEAR_SOLS:.0f} sols)")
    yticks = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    yticks = [t for t in yticks if t <= period_sols.max()]
    ax2.set_yticks(np.log2(yticks))
    ax2.set_yticklabels(yticks)
    ax2.set_ylabel("Period (sols)")
    ax2.set_xlabel("SOL (mission day)")
    ax2.set_title(f"Wavelet power spectrum, detector {detector}, sol {sol_axis[0]}-{sol_axis[-1]} "
                  f"({n_cycles:.1f} Martian years) -- Morlet, white=95% signif, hatched=COI")
    ax2.legend(loc="upper right")
    fig2.colorbar(cf, ax=ax2, label="Wavelet power (normalized units)")
    fig2.tight_layout()
    out_png2 = HERE / f"wavelet_{detector}_extended_full_history.png"
    fig2.savefig(out_png2, dpi=200, bbox_inches="tight")
    plt.close(fig2)
    print(f"  saved {out_png2}")

    # check power specifically in the annual band vs surrounding bands
    annual_mask = (period_sols > MARS_YEAR_SOLS * 0.85) & (period_sols < MARS_YEAR_SOLS * 1.15)
    annual_power_frac_significant = (sig95[annual_mask, :] > 1).mean()
    print(f"  fraction of annual-band power exceeding 95% significance: {annual_power_frac_significant*100:.2f}%")

    return sol_axis, y, trend, seasonal, resid


if __name__ == "__main__":
    df = load_master()
    for det in ["B", "E"]:
        run_for_detector(df, det)
