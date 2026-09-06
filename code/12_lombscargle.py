"""
12_lombscargle.py

Third, independent test for a seasonal (annual, ~668.6-sol) periodicity in the
RAD dose-rate record, complementing the classical decomposition (phase-averaged,
sensitive to the exact period assumed) and the wavelet transform (localized in
time, assessed against a red-noise null). The Lomb-Scargle periodogram tests
the full record globally, with a formal false-alarm-probability significance
level, independent of any phase-binning assumption -- closes the residual
"maybe the phase-averaging method missed it" objection.

Run on the SOLAR-REGRESSION RESIDUAL (dose minus the NM-based fit from
11_nm_lag_hysteresis.py), i.e. after removing the dominant secular/solar-cycle
trend, since a periodogram on the raw series would be swamped by that trend's
power at low frequencies.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from astropy.timeseries import LombScargle
from pathlib import Path

HERE = Path(__file__).parent
MASTER = HERE / "mslrad_master_sol_series.csv"
NM_CSV = HERE / "nm_data" / "oulu_daily.csv"
MARS_YEAR_SOLS = 668.6
SMOOTH_WINDOW_SOLS = 387


def load_master():
    return pd.read_csv(MASTER, parse_dates=["date"])


def load_nm_on_sol_grid(dates):
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
    return pd.Series(np.where(mask, np.nan, y)).interpolate(limit_direction="both").values


def run_detector(df, detector, nm_grid):
    print(f"\n{'='*70}\nDetector {detector}\n{'='*70}")
    dose_raw = df[f"dose_{detector}_mean"].values
    dose, first = interpolate_gaps(dose_raw)
    n = len(dose)
    sol_axis = df["sol"].values[first:first + n].astype(float)
    dose_clean = clean_spe(dose)

    nm = nm_grid[first:first + n]
    nm_smooth = pd.Series(nm).rolling(SMOOTH_WINDOW_SOLS, center=True, min_periods=1).mean().values
    X = sm.add_constant(nm_smooth)
    model = sm.OLS(dose_clean, X).fit()
    resid = dose_clean - model.predict(X)
    print(f"solar-trend (NM-based) fit: R2={model.rsquared:.3f}")

    # periods to scan: a few sols up to ~1500 sols (over 2 Martian years), log-spaced
    periods = np.geomspace(4, 1500, 2000)
    freqs = 1.0 / periods
    ls = LombScargle(sol_axis, resid)
    power = ls.power(freqs)

    # false alarm probability levels -- for the GLOBAL scan (2000 frequencies, few sols to
    # 1500 sols), correcting for the number of independent frequencies tested (Baluev 2008)
    fap_levels = [0.1, 0.05, 0.01]
    fap_power = ls.false_alarm_level(fap_levels, method="baluev")

    # power at the SINGLE, pre-specified annual frequency (1/668.6 sols) -- this is a targeted
    # test of one externally motivated hypothesis (Guo et al. 2015's seasonal-pressure claim),
    # not a peak selected by scanning, so no multiple-testing correction is applied: the
    # classical single-frequency false alarm probability for astropy's "standard"-normalized
    # Lomb-Scargle power is exp(-power) (Scargle 1982).
    annual_power = ls.power(np.array([1.0 / MARS_YEAR_SOLS]))[0]
    annual_fap_single = np.exp(-annual_power)
    print(f"power at the single annual frequency (1/{MARS_YEAR_SOLS:.1f} sols): {annual_power:.4f}, "
          f"single-frequency FAP={annual_fap_single:.3f}")

    # global strongest peak anywhere in the scanned range -- this DOES need the multi-frequency
    # correction, since it was found by scanning (report for context, not as the primary test)
    i_max = np.argmax(power)
    print(f"[context only, not the annual test] strongest peak anywhere in scan: "
          f"period={periods[i_max]:.1f} sols, power={power[i_max]:.4f}, "
          f"scan-corrected FAP={ls.false_alarm_probability(power[i_max], method='baluev'):.3e}")
    print(f"scan-corrected FAP thresholds (power level): 10%={fap_power[0]:.4f}, 5%={fap_power[1]:.4f}, 1%={fap_power[2]:.4f}")

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.semilogx(periods, power, color="navy", lw=0.8)
    for lvl, p in zip(fap_levels, fap_power):
        ax.axhline(p, ls="--", lw=0.8, color="gray")
        ax.text(periods[-1], p, f" FAP={lvl:.0%}", va="bottom", ha="right", fontsize=7, color="gray")
    ax.axvline(MARS_YEAR_SOLS, color="red", ls=":", lw=1.2, label=f"1 Martian year ({MARS_YEAR_SOLS:.0f} sols)")
    ax.annotate("rises toward the longest periods scanned:\nresidual trend curvature (a single OLS line\ncan't fully capture the trend's real shape),\nnot a periodic signal",
                xy=(1300, power.max() * 0.85), fontsize=7, color="dimgray", ha="right")
    ax.set_xlabel("Period (sols)")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title(f"Detector {detector}: Lomb-Scargle periodogram of the solar-regression residual")
    ax.legend(loc="upper right")
    fig.tight_layout()
    out_png = HERE / f"lombscargle_{detector}.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    return {
        "detector": detector, "annual_power": annual_power, "annual_fap_single": annual_fap_single,
        "global_peak_period": periods[i_max], "global_peak_power": power[i_max],
        "global_peak_fap_scan_corrected": ls.false_alarm_probability(power[i_max], method="baluev"),
    }


if __name__ == "__main__":
    df = load_master()
    nm_grid = load_nm_on_sol_grid(df["date"])
    results = [run_detector(df, det, nm_grid) for det in ["B", "E"]]
    out = pd.DataFrame(results)
    out.to_csv(HERE / "lombscargle_summary.csv", index=False)
    print("\n" + out.to_string(index=False))
