"""
10_solar_pressure_regression.py

Fecha o ponto (1) do plano de fortalecimento do paper: em vez de so citar
Harri et al. (2014) em prosa e testar sazonalidade "as cegas" na dose bruta,
regride a dose diretamente contra (a) atividade solar continua (numero de
manchas solares suavizado, SILSO observado ate 2026-08 -- nao mais um corte
binario Period1/Period2) e (b) pressao real do REMS (extraida em
06_extract_rems_pressure.py), testando o RESIDUO por periodicidade sazonal.

Controle positivo: a propria pressao passa pela mesma decomposicao classica
usada na dose. A pressao tem um ciclo sazonal real e bem documentado (CO2
condensando/sublimando nas calotas polares) -- se o metodo achar uma razao
sazonal/residuo grande na pressao mas pequena no residuo da dose, isso
demonstra que a ausencia de sazonalidade na dose NAO e um artefato de um
metodo insensivel, e sim um resultado negativo real.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pathlib import Path

HERE = Path(__file__).parent
MASTER = HERE / "mslrad_master_sol_series.csv"
MARS_YEAR_SOLS = 668.6
SMOOTH_WINDOW_SOLS = 387  # ~13 Earth months in sols, the standard "smoothed SSN" window

# SOL 4190 (~2024-05-14) coincides with the well-documented May 2024 solar storm
# (AR3664, X8.7 flare on 2024-05-14) -- NASA/JPL reported it as the largest
# radiation surge RAD had recorded since landing (2012). Per-sol mean dose here
# is B=71.7, E=78.4 uGy/hour vs a record mean of ~9.8-9.9 -- an SPE transient,
# not part of the smooth secular/solar-cycle trend this section tests. Following
# the same convention already used in Table 1 (windows flagged "2014 SPE",
# "2017 SPE" kept in descriptive stats but treated separately from the trend
# discussion), it is excluded from the OLS solar-regression fit below; the
# fit including it is also reported, for transparency.
SPE_EXCLUDED_SOLS = [4190]


def load_master():
    return pd.read_csv(MASTER, parse_dates=["date"])


def interpolate_gaps(y: np.ndarray):
    s = pd.Series(y)
    first, last = s.first_valid_index(), s.last_valid_index()
    s = s.iloc[first:last + 1]
    frac_missing = s.isna().mean()
    return s.interpolate(method="linear", limit_direction="both").values, frac_missing, first


def flag_spe_outlier_sols(y: np.ndarray, window: int = 21, mad_threshold: float = 8.0) -> np.ndarray:
    """Same robust-outlier flag as 09_extended_analysis.py: sols departing from
    a short rolling median by more than mad_threshold robust MADs -- needed
    here too, since seasonal_to_residual_ratio() below runs the same
    phase-averaged classical_decompose(), which is not robust to SPE spikes
    (one outlier smears into a fake repeated "seasonal" spike; see that
    script's docstring). Excluding only the single largest SPE (sol 4190, in
    the OLS fit above) is not enough for THIS ratio -- the smaller, still-present
    SPEs (~8-9 more, matching Table 1's flagged windows plus new solar-cycle-25
    events) would otherwise still inflate it."""
    s = pd.Series(y)
    roll_med = s.rolling(window, center=True, min_periods=1).median()
    dev = (s - roll_med).abs()
    mad = 1.4826 * dev.median()
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


def seasonal_to_residual_ratio(y, period=round(MARS_YEAR_SOLS)):
    y_i, frac_missing, first = interpolate_gaps(y)
    spe_mask = flag_spe_outlier_sols(y_i)
    y_clean = np.where(spe_mask, np.nan, y_i)
    y_clean = pd.Series(y_clean).interpolate(limit_direction="both").values
    _, seasonal, resid = classical_decompose(y_clean, period)
    ratio = np.nanstd(seasonal) / np.nanstd(resid)
    return ratio, seasonal, resid, y_i, first, frac_missing


def cross_correlate(a, b, max_lag):
    """Pearson correlation of a[t] vs b[t+lag] for lag in [-max_lag, max_lag],
    over the overlapping region, both z-scored. Returns (lags, r)."""
    a = (a - np.nanmean(a)) / np.nanstd(a)
    b = (b - np.nanmean(b)) / np.nanstd(b)
    lags = np.arange(-max_lag, max_lag + 1)
    rs = []
    for lag in lags:
        if lag >= 0:
            aa, bb = a[lag:], b[:len(b) - lag]
        else:
            aa, bb = a[:len(a) + lag], b[-lag:]
        mask = ~np.isnan(aa) & ~np.isnan(bb)
        if mask.sum() < 30:
            rs.append(np.nan)
            continue
        rs.append(np.corrcoef(aa[mask], bb[mask])[0, 1])
    return lags, np.array(rs)


def run_detector(df, detector):
    print(f"\n{'='*70}\nDetector {detector}\n{'='*70}")
    dose_raw = df[f"dose_{detector}_mean"].values
    dose, frac_missing, first = interpolate_gaps(dose_raw)
    n = len(dose)
    sol_axis = df["sol"].values[first:first + n]

    ssn = df["ssn"].values[first:first + n]
    ssn_smooth = pd.Series(ssn).rolling(SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values

    # --- (a) regress dose against smoothed solar activity, continuous, not a binary split ---
    X = sm.add_constant(ssn_smooth)
    model_full = sm.OLS(dose, X).fit()
    print(f"OLS: dose_{detector} ~ smoothed_SSN (all sols)")
    print(f"  slope={model_full.params[1]:.5f} (uGy/h per SSN unit), p={model_full.pvalues[1]:.3e}, "
          f"R2={model_full.rsquared:.4f}")

    spe_mask = ~np.isin(sol_axis, SPE_EXCLUDED_SOLS)
    model = sm.OLS(dose[spe_mask], X[spe_mask]).fit()
    print(f"  excluding SPE sol(s) {SPE_EXCLUDED_SOLS}: slope={model.params[1]:.5f}, "
          f"p={model.pvalues[1]:.3e}, R2={model.rsquared:.4f}  <- reported as primary fit")

    resid_full = dose - model_full.predict(X)
    resid = np.where(spe_mask, dose - model.predict(X), np.nan)

    # --- positive control: does the SAME method find the known real seasonal cycle in pressure? ---
    pressure_raw = df["pressure_pa_mean"].values
    pressure_aligned = pressure_raw[first:first + n]
    p_ratio, p_seasonal, p_resid, p_interp, p_first_rel, p_frac_missing = seasonal_to_residual_ratio(pressure_aligned)
    print(f"  [positive control] pressure seasonal/residual ratio = {p_ratio:.3f} "
          f"({p_frac_missing*100:.1f}% of pressure series interpolated)")

    # --- (b) test the solar-regression RESIDUAL for a seasonal (annual) component ---
    r_ratio, r_seasonal, r_resid2, r_interp, r_first_rel, r_frac_missing = seasonal_to_residual_ratio(resid)
    print(f"  dose residual (after removing solar trend) seasonal/residual ratio = {r_ratio:.3f}")
    print(f"  => pressure shows {p_ratio/r_ratio:.1f}x more seasonal structure than the dose residual "
          f"(same method, same time axis)")

    # --- (c) cross-correlate dose residual directly against real REMS pressure ---
    lags, rs = cross_correlate(resid, pressure_aligned, max_lag=400)
    best_i = np.nanargmax(np.abs(rs))
    print(f"  cross-correlation (dose residual vs REMS pressure, lags -400..+400 sols): "
          f"max |r|={abs(rs[best_i]):.4f} at lag={lags[best_i]} sols")

    fig, axes = plt.subplots(3, 1, figsize=(13, 8))
    axes[0].plot(sol_axis, dose, lw=0.4, color="navy", label="observed dose")
    axes[0].plot(sol_axis, model.predict(X), lw=1.2, color="red",
                 label=f"OLS fit vs. smoothed SSN (R²={model.rsquared:.2f}, SPE sol excluded)")
    for spe_sol in SPE_EXCLUDED_SOLS:
        if sol_axis.min() <= spe_sol <= sol_axis.max():
            axes[0].annotate(f"SOL {spe_sol}\n(May 2024 SPE)", xy=(spe_sol, dose[sol_axis == spe_sol][0]),
                              xytext=(spe_sol - 500, dose.max() * 0.7), fontsize=7, color="dimgray",
                              arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.7))
    axes[0].set_ylabel("Dose\n(uGy/hour)")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].set_title(f"Detector {detector}: dose vs. smoothed solar activity, and residual seasonal test")

    axes[1].plot(sol_axis, resid, lw=0.3, color="gray", label="dose residual (obs - solar fit)")
    axes[1].set_ylabel("Residual\n(uGy/hour)")
    axes[1].legend(loc="upper right", fontsize=8)

    ax2b = axes[1].twinx()
    ax2b.plot(sol_axis, pressure_aligned, lw=0.5, color="teal", alpha=0.6, label="REMS pressure")
    ax2b.set_ylabel("Pressure (Pa)", color="teal")

    axes[2].plot(lags, rs, color="darkorange")
    axes[2].axvline(0, color="k", lw=0.5, ls=":")
    axes[2].axhline(0, color="k", lw=0.5)
    axes[2].set_xlabel("Lag (sols); positive = pressure leads dose residual")
    axes[2].set_ylabel("Cross-correlation r")
    fig.tight_layout()
    out_png = HERE / f"solar_pressure_regression_{detector}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    return {
        "detector": detector, "slope": model.params[1], "p_value": model.pvalues[1],
        "r2": model.rsquared, "r2_all_sols_incl_spe": model_full.rsquared,
        "pressure_seasonal_ratio": p_ratio, "residual_seasonal_ratio": r_ratio,
        "max_abs_xcorr": abs(rs[best_i]), "xcorr_lag": lags[best_i],
    }


if __name__ == "__main__":
    df = load_master()
    results = [run_detector(df, det) for det in ["B", "E"]]
    out = pd.DataFrame(results)
    out.to_csv(HERE / "solar_pressure_regression_summary.csv", index=False)
    print("\n" + out.to_string(index=False))
