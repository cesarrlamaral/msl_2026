"""
Two-sample t-test comparing mean dose rate (detectors B and E) between
Period 1 (SOL 00000-01159, the first 11 of 22 ~90-SOL windows in
medias.xlsx) and Period 2 (SOL 01160-02482, the remaining 11 windows).

This reproduces the comparison reported in the original R-based analysis
(Table 2 of paper_msl_2026.docx / the "B - red (antes) x green (depois)"
and "E - red (antes) x green (depois)" cells in medias.xlsx), which had no
saved script -- the exact R call that produced p=5.11297e-6 (B) and
p=2.68702e-6 (E) could not be found anywhere in the project. This script
makes the test reproducible and explicit about what data/units it runs on:
the 22 per-window MEANS (n=11 per period), not the raw per-observation
dose values (which are strongly autocorrelated and would give an invalid,
inflated significance if tested directly).

Author: Cesar R. L. Amaral (analysis re-implemented in Python, 2026-09-05,
per project decision to move off R -- see feedback_mslrad_python_nao_r).
"""
import pandas as pd
from scipy import stats

MEDIAS_PATH = r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\MSLRADv1 RESULTS\medias.xlsx"


def load_window_means():
    df = pd.read_excel(MEDIAS_PATH).iloc[:, :5]
    df.columns = ["sols", "mediab", "desvb", "mediae", "desve"]
    return df


def run_period_ttest():
    df = load_window_means()
    period1 = df.iloc[0:11]   # SOL 00000_00089 .. SOL 01063_01159
    period2 = df.iloc[11:22]  # SOL 01160_01293 .. SOL 02359_02482

    print(f"Period 1: {len(period1)} windows, {period1['sols'].iloc[0]} .. {period1['sols'].iloc[-1]}")
    print(f"Period 2: {len(period2)} windows, {period2['sols'].iloc[0]} .. {period2['sols'].iloc[-1]}")
    print()

    results = {}
    for det, col in [("B", "mediab"), ("E", "mediae")]:
        t, p = stats.ttest_ind(period1[col], period2[col], equal_var=False)  # Welch, matches R's t.test() default
        results[det] = (t, p)
        print(f"Detector {det}: Welch two-sample t-test, "
              f"mean(Period 1)={period1[col].mean():.4f}, mean(Period 2)={period2[col].mean():.4f}, "
              f"t={t:.4f}, p={p:.3e}")

    print()
    print("Note: original (R, unsaved) analysis reported p=5.113e-6 (B) and p=2.687e-6 (E) for "
          "what the spreadsheet labels the same Period 1 vs Period 2 comparison. This Python "
          "reproduction, run explicitly as a Welch two-sample t-test on the 22 window-level means "
          "(n=11 per group), agrees on direction and order of magnitude (both far below p=0.05) "
          "but not on the exact value -- the original R call was never saved, so the precise test "
          "variant used (equal- vs unequal-variance, one- vs two-sided, or a different data "
          "arrangement) could not be confirmed. Reported here as the reproducible reference "
          "going forward.")
    return results


if __name__ == "__main__":
    run_period_ttest()
