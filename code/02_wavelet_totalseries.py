"""
Continuous wavelet transform (Morlet, Torrence & Compo 1998) of the MSLRAD
TOTALSERIES dose-rate record (detector E, SOL 00000-02358), replacing the R
WaveletComp pass that used upperPeriod=500 observations (~5 sols) -- far too
short to reveal a Martian seasonal (~668.6 sol) or multi-year (solar-cycle
adjacent) signal in the combined 7-year series. Here the scale range is
extended to span from a few sols up to roughly 1/3 of the record length
(the largest period for which the cone of influence still covers a
meaningful central portion of the data).

Author: Cesar R. L. Amaral (analysis re-implemented in Python, 2026-09-05,
per project decision to move off R -- see feedback_mslrad_python_nao_r).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pycwt as wavelet

DATA_DIR = r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\MSLRADv1 TOTALSERIES RESULTS"
OUT_DIR = r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\python_reanalysis"

SOL_SPAN = 2358
MARS_YEAR_SOLS = 668.6


def load_series(detector):
    path = f"{DATA_DIR}\\MSLRAD_TOTALSERIES_{detector}.txt"
    s = pd.read_csv(path, header=0).iloc[:, 0].values.astype(float)
    return s


def run_wavelet(detector):
    y = load_series(detector)
    n = len(y)
    obs_per_sol = n / SOL_SPAN
    dt = 1.0  # 1 time unit = 1 observation
    sol_axis = np.arange(n) / obs_per_sol

    # Standardize (remove mean, normalize by std) -- standard practice for CWT
    y0 = y - y.mean()
    std = y0.std()
    y_norm = y0 / std

    mother = wavelet.Morlet(6)
    s0 = 2 * dt  # smallest resolvable scale (~2 observations)
    dj = 1 / 12  # 12 sub-octaves per octave (finer than R script's 1/50 was for a much
                 # narrower range; 1/12 keeps output size manageable over this wide range)
    # largest scale: allow up to ~1/3 of the record so the cone of influence
    # still covers the central ~1/3 of the series at the longest periods probed
    J = int(np.log2(n / 3 / s0) / dj)

    print(f"[{detector}] n={n}, obs/sol={obs_per_sol:.3f}, scanning up to "
          f"~{(n/3)/obs_per_sol:.0f} sols period, {J+1} scales")

    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(y_norm, dt, dj, s0, J, mother)
    power = (np.abs(wave)) ** 2
    period = 1 / freqs

    # significance test against a red-noise background (standard Torrence & Compo)
    alpha, _, _ = wavelet.ar1(y_norm)
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                              significance_level=0.95,
                                              wavelet=mother)
    sig95 = np.ones([1, n]) * signif[:, None]
    sig95 = power / sig95

    period_sols = period / obs_per_sol
    coi_sols = coi / obs_per_sol

    fig, ax = plt.subplots(figsize=(13, 6))
    levels = np.linspace(0, np.percentile(power, 99.5), 30)
    cf = ax.contourf(sol_axis, np.log2(period_sols), np.log2(power) if False else power,
                      levels=levels, extend="both", cmap="viridis")
    ax.contour(sol_axis, np.log2(period_sols), sig95, [1], colors="white", linewidths=1.2)
    ax.fill(np.concatenate([sol_axis, sol_axis[-1:], sol_axis[:1]]),
            np.concatenate([np.log2(coi_sols), [np.log2(period_sols.max())], [np.log2(period_sols.max())]]),
            "k", alpha=0.4, hatch="/")
    ax.axhline(np.log2(MARS_YEAR_SOLS), color="red", ls="--", lw=1,
               label=f"1 Martian year ({MARS_YEAR_SOLS:.0f} sols)")
    yticks = [8, 16, 32, 64, 128, 256, 512, 1024]
    yticks = [t for t in yticks if t <= period_sols.max()]
    ax.set_yticks(np.log2(yticks))
    ax.set_yticklabels(yticks)
    ax.set_ylabel("Period (sols)")
    ax.set_xlabel("SOL (mission day)")
    ax.set_title(f"Wavelet power spectrum, detector {detector}, TOTALSERIES "
                 f"(SOL 0-{SOL_SPAN}) -- Morlet, white contour = 95% significance "
                 f"vs. red noise, hatched = cone of influence")
    ax.legend(loc="upper right")
    fig.colorbar(cf, ax=ax, label="Wavelet power (normalized units)")
    fig.tight_layout()

    out_png = f"{OUT_DIR}\\Wavelet_TOTALSERIES_{detector}_python_fullrange.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_png}")

    return period_sols, power, sig95


if __name__ == "__main__":
    run_wavelet("E")
    run_wavelet("B")
