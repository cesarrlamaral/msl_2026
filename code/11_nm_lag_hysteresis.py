"""
11_nm_lag_hysteresis.py

Addresses three of the "high-impact" gaps identified in critical review of the
extended paper (2026-09-06):

(1) Autocorrelation-robust significance: the Section 3.4 OLS regression
    (dose ~ smoothed SSN) reports p-values as if n=4843 independent sols, but
    both series are heavily autocorrelated (dose by physical persistence, SSN
    by construction via a 387-sol smoothing window) -- the true effective N is
    far smaller. Refit here with Newey-West (HAC) standard errors.

(2) Lag/hysteresis test: the SSN regression used an *instantaneous* smoothed
    SSN, never scanning for a lag -- unlike the pressure cross-correlation in
    10_solar_pressure_regression.py, which did scan lags. GCR modulation is
    known (from Earth neutron-monitor literature) to lag solar activity, often
    asymmetrically between the rising and declining phases of a solar cycle
    (a "hysteresis loop" in GCR-vs-activity plots). This record now spans a
    FULL cycle reversal (24 declining -> minimum -> 25 rising) for the first
    time in the MSL/RAD literature -- no prior published RAD paper (Guo 2015,
    Guo 2021, Ehresmann 2023) had more than a monotonic decline, so this test
    was not previously possible.

(3) Physical GCR proxy: replaces/supplements the sunspot-number proxy with the
    Oulu neutron monitor daily count rate (NMDB, corr_for_efficiency),
    2012-08-01 to 2026-08-31 -- a direct, physically motivated measurement of
    the same heliospheric GCR modulation Mars-surface dose responds to, and
    enables a real quantitative Earth-vs-Mars modulation-amplitude comparison.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pathlib import Path

HERE = Path(__file__).parent
MASTER = HERE / "mslrad_master_sol_series.csv"
NM_CSV = HERE / "nm_data" / "oulu_daily.csv"
LANDING_DATE = pd.Timestamp("2012-08-06")
SOL_TO_DAYS = 1.0274912
SMOOTH_WINDOW_SOLS = 387
LAG_SMOOTH_WINDOW_SOLS = 27  # ~Bartels solar-rotation period, standard for lag/cross-correlation studies
                              # (387 sols is too wide to resolve a realistic month-scale lag -- it would
                              # smear/alias any lag shorter than the smoothing window itself)
HAC_MAXLAGS = 400  # matches the smoothing window; conservative given ~400-800 sol autocorrelation scale
SPE_SOLS = [4190]
TREND_TURNING_SOL = 2600  # approx cycle24-decline / cycle25-rise inflection, from Figure 5 (decomposition trend)


def load_master():
    return pd.read_csv(MASTER, parse_dates=["date"])


def load_nm_on_sol_grid(dates: pd.Series) -> np.ndarray:
    nm = pd.read_csv(NM_CSV, parse_dates=["date"])
    x = nm["date"].values.astype("datetime64[ns]").astype("int64")
    y = nm["count_rate"].values
    xq = dates.values.astype("datetime64[ns]").astype("int64")
    return np.interp(xq, x, y)


def interpolate_gaps(y):
    s = pd.Series(y)
    first, last = s.first_valid_index(), s.last_valid_index()
    s = s.iloc[first:last + 1]
    return s.interpolate(limit_direction="both").values, first


def flag_spe_outlier_sols(y, window=21, mad_threshold=8.0):
    s = pd.Series(y)
    roll_med = s.rolling(window, center=True, min_periods=1).median()
    dev = (s - roll_med).abs()
    mad = 1.4826 * dev.median()
    if mad == 0 or np.isnan(mad):
        return np.zeros(len(y), dtype=bool)
    return (dev > mad_threshold * mad).values


def clean_spe(y):
    mask = flag_spe_outlier_sols(y)
    y2 = np.where(mask, np.nan, y)
    return pd.Series(y2).interpolate(limit_direction="both").values


def lag_scan(dose, regressor, sol_axis, max_lag=800, step=5):
    """Slide regressor relative to dose; positive lag = regressor LEADS dose
    (i.e. dose responds to regressor lag sols earlier). Returns lags, r2 of
    simple OLS fit at each lag."""
    lags = np.arange(-max_lag, max_lag + 1, step)
    r2s = []
    for lag in lags:
        if lag >= 0:
            d = dose[lag:]
            r = regressor[:len(regressor) - lag] if lag > 0 else regressor
        else:
            d = dose[:len(dose) + lag]
            r = regressor[-lag:]
        n = min(len(d), len(r))
        d, r = d[:n], r[:n]
        if n < 200:
            r2s.append(np.nan)
            continue
        X = sm.add_constant(r)
        m = sm.OLS(d, X).fit()
        r2s.append(m.rsquared)
    return lags, np.array(r2s)


def run_detector(df, detector, nm_on_grid):
    print(f"\n{'='*70}\nDetector {detector}\n{'='*70}")
    dose_raw = df[f"dose_{detector}_mean"].values
    dose, first = interpolate_gaps(dose_raw)
    n = len(dose)
    sol_axis = df["sol"].values[first:first + n]
    dose_clean = clean_spe(dose)

    ssn = df["ssn"].values[first:first + n]
    ssn_smooth = pd.Series(ssn).rolling(SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values
    nm = nm_on_grid[first:first + n]
    nm_smooth = pd.Series(nm).rolling(SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values
    nm_lightsmooth = pd.Series(nm).rolling(LAG_SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values
    dose_lightsmooth = pd.Series(dose_clean).rolling(LAG_SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values

    # --- (1) HAC-corrected refit of the original SSN regression ---
    X_ssn = sm.add_constant(ssn_smooth)
    m_ssn_ols = sm.OLS(dose_clean, X_ssn).fit()
    m_ssn_hac = sm.OLS(dose_clean, X_ssn).fit(cov_type="HAC", cov_kwds={"maxlags": HAC_MAXLAGS})
    print(f"[SSN, instantaneous] OLS: slope={m_ssn_ols.params[1]:.5f} p={m_ssn_ols.pvalues[1]:.3e} R2={m_ssn_ols.rsquared:.3f}")
    print(f"[SSN, instantaneous] HAC (maxlags={HAC_MAXLAGS}): p={m_ssn_hac.pvalues[1]:.3e}, "
          f"SE ratio (HAC/OLS)={m_ssn_hac.bse[1]/m_ssn_ols.bse[1]:.1f}x")

    # --- (2) lag scan: dose vs NM count rate (physical proxy), lightly smoothed
    # (27-sol / Bartels-rotation scale) so a realistic month-scale lag is not
    # washed out by the 387-sol window used for the headline regression ---
    lags, r2s = lag_scan(dose_lightsmooth, nm_lightsmooth, sol_axis, max_lag=800, step=5)
    best_i = np.nanargmax(r2s)
    best_lag = lags[best_i]
    print(f"[NM lag scan, {LAG_SMOOTH_WINDOW_SOLS}-sol smoothing] best lag = {best_lag} sols (R2={r2s[best_i]:.3f}), "
          f"vs instantaneous (lag=0) R2={r2s[np.where(lags==0)[0][0]]:.3f}")

    # refit at best lag with HAC
    if best_lag >= 0:
        d_fit, r_fit, sol_fit = dose_clean[best_lag:], nm_smooth[:len(nm_smooth)-best_lag], sol_axis[best_lag:]
    else:
        d_fit, r_fit, sol_fit = dose_clean[:len(dose_clean)+best_lag], nm_smooth[-best_lag:], sol_axis[:len(sol_axis)+best_lag]
    nfit = min(len(d_fit), len(r_fit))
    d_fit, r_fit, sol_fit = d_fit[:nfit], r_fit[:nfit], sol_fit[:nfit]
    X_nm = sm.add_constant(r_fit)
    m_nm_ols = sm.OLS(d_fit, X_nm).fit()
    m_nm_hac = sm.OLS(d_fit, X_nm).fit(cov_type="HAC", cov_kwds={"maxlags": HAC_MAXLAGS})
    print(f"[NM, lag={best_lag}] OLS: slope={m_nm_ols.params[1]:.5f} R2={m_nm_ols.rsquared:.3f}; "
          f"HAC p={m_nm_hac.pvalues[1]:.3e}")

    # instantaneous (lag=0) NM fit too, for direct comparison to SSN
    X_nm0 = sm.add_constant(nm_smooth)
    m_nm0_ols = sm.OLS(dose_clean, X_nm0).fit()
    m_nm0_hac = sm.OLS(dose_clean, X_nm0).fit(cov_type="HAC", cov_kwds={"maxlags": HAC_MAXLAGS})
    print(f"[NM, lag=0] OLS: slope={m_nm0_ols.params[1]:.5f} R2={m_nm0_ols.rsquared:.3f}; HAC p={m_nm0_hac.pvalues[1]:.3e}")

    # --- (3) hysteresis: split at trend turning point, compare phase-specific fits ---
    decl_mask = sol_fit < TREND_TURNING_SOL
    rise_mask = ~decl_mask
    m_decl = sm.OLS(d_fit[decl_mask], X_nm[decl_mask]).fit()
    m_rise = sm.OLS(d_fit[rise_mask], X_nm[rise_mask]).fit()
    print(f"[hysteresis @ best lag] declining phase (sol<{TREND_TURNING_SOL}): "
          f"slope={m_decl.params[1]:.5f}, intercept={m_decl.params[0]:.3f}, R2={m_decl.rsquared:.3f}, n={decl_mask.sum()}")
    print(f"[hysteresis @ best lag] rising phase (sol>={TREND_TURNING_SOL}): "
          f"slope={m_rise.params[1]:.5f}, intercept={m_rise.params[0]:.3f}, R2={m_rise.rsquared:.3f}, n={rise_mask.sum()}")

    # formal test of slope difference: dose ~ NM + phase + NM:phase, HAC-robust SE on the interaction term
    phase = rise_mask.astype(float)
    X_inter = np.column_stack([np.ones(nfit), r_fit, phase, r_fit * phase])
    m_inter = sm.OLS(d_fit, X_inter).fit(cov_type="HAC", cov_kwds={"maxlags": HAC_MAXLAGS})
    slope_diff, p_diff = m_inter.params[3], m_inter.pvalues[3]
    print(f"[hysteresis significance] NM x phase interaction (HAC): "
          f"slope_diff={slope_diff:.5f} (rising - declining), p={p_diff:.3e}")

    # --- plot: lag scan curve + hysteresis loop ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(lags, r2s, color="navy")
    axes[0].axvline(best_lag, color="red", ls="--", lw=1, label=f"best lag = {best_lag} sols")
    axes[0].axvline(0, color="gray", lw=0.5)
    axes[0].set_xlabel("Lag (sols); positive = NM count rate leads dose")
    axes[0].set_ylabel("R² (dose ~ NM count rate)")
    axes[0].set_title(f"Detector {detector}: lag scan, dose vs. Oulu NM count rate")
    axes[0].legend()

    sc = axes[1].scatter(r_fit, d_fit, c=sol_fit, cmap="coolwarm", s=4, alpha=0.5)
    axes[1].plot(r_fit[decl_mask][np.argsort(r_fit[decl_mask])],
                 m_decl.predict(X_nm[decl_mask])[np.argsort(r_fit[decl_mask])], color="blue", lw=2,
                 label=f"declining phase (sol<{TREND_TURNING_SOL})")
    axes[1].plot(r_fit[rise_mask][np.argsort(r_fit[rise_mask])],
                 m_rise.predict(X_nm[rise_mask])[np.argsort(r_fit[rise_mask])], color="red", lw=2,
                 label=f"rising phase (sol≥{TREND_TURNING_SOL})")
    axes[1].set_xlabel(f"Oulu NM count rate (%, lag={best_lag} sols)")
    axes[1].set_ylabel("Dose rate (uGy/hour)")
    axes[1].set_title("Hysteresis: dose vs. GCR proxy, colored by sol")
    axes[1].legend(fontsize=8)
    fig.colorbar(sc, ax=axes[1], label="SOL")
    fig.tight_layout()
    out_png = HERE / f"lag_hysteresis_{detector}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    return {
        "detector": detector,
        "ssn_slope": m_ssn_ols.params[1], "ssn_r2": m_ssn_ols.rsquared,
        "ssn_p_ols": m_ssn_ols.pvalues[1], "ssn_p_hac": m_ssn_hac.pvalues[1],
        "nm_lag0_slope": m_nm0_ols.params[1], "nm_lag0_r2": m_nm0_ols.rsquared, "nm_lag0_p_hac": m_nm0_hac.pvalues[1],
        "nm_best_lag": best_lag, "nm_best_r2": m_nm_ols.rsquared, "nm_best_p_hac": m_nm_hac.pvalues[1],
        "decl_slope": m_decl.params[1], "decl_r2": m_decl.rsquared,
        "rise_slope": m_rise.params[1], "rise_r2": m_rise.rsquared,
        "slope_diff": slope_diff, "slope_diff_p_hac": p_diff,
        "lag_scan_smoothing_sols": LAG_SMOOTH_WINDOW_SOLS,
    }


if __name__ == "__main__":
    df = load_master()
    nm_grid = load_nm_on_sol_grid(df["date"])
    results = [run_detector(df, det, nm_grid) for det in ["B", "E"]]
    out = pd.DataFrame(results)
    out.to_csv(HERE / "nm_lag_hysteresis_summary.csv", index=False)
    print("\n" + out.to_string(index=False))
