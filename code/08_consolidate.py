"""
08_consolidate.py

Consolida as extracoes byte-exatas (RAD dose + REMS pressao, sols 1-4843) em
uma unica serie por sol, com data real e numero de manchas solares (SILSO,
observado, atualizado ate 2026-08) -- a base para a reanalise estatistica
completa (decomposicao/wavelet na serie estendida + regressao solar/pressao).

Saida: mslrad_master_sol_series.csv, colunas:
  sol, date, dose_B_mean, dose_B_n, dose_E_mean, dose_E_n,
  pressure_pa_mean, pressure_n, ssn
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
LANDING_DATE = pd.Timestamp("2012-08-06")  # MSL landing, SOL 0
SOL_TO_DAYS = 1.0274912  # 1 Martian sol in Earth days

SSN_PATH = HERE.parent / "MSLRAD BIBLIOS" / "paper msl rad" / "SN_m_tot_V2.0_UPDATED.txt"


def load_dose() -> pd.DataFrame:
    files = sorted(glob.glob(str(HERE / "extended_series" / "*.csv")))
    print(f"loading {len(files)} dose CSVs")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.drop_duplicates(subset=["sol", "obs_index", "detector"])
    agg = df.pivot_table(index="sol", columns="detector", values="value", aggfunc=["mean", "count"])
    agg.columns = [f"dose_{det}_{stat}" for stat, det in agg.columns]
    agg = agg.rename(columns={"dose_B_count": "dose_B_n", "dose_E_count": "dose_E_n"})
    return agg.reset_index()


def load_pressure() -> pd.DataFrame:
    files = sorted(glob.glob(str(HERE / "pressure_series" / "*.csv")))
    print(f"loading {len(files)} pressure CSVs")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.dropna(subset=["pressure_pa"])
    g = df.groupby("sol")["pressure_pa"].agg(["mean", "count"]).reset_index()
    g.columns = ["sol", "pressure_pa_mean", "pressure_n"]
    return g


def load_ssn_monthly() -> pd.DataFrame:
    # SILSO format: year month decimal_year SSN SSN_std n_obs [marker]
    cols = ["year", "month", "dec_year", "ssn", "ssn_std", "nobs", "flag"]
    rows = []
    with open(SSN_PATH, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            rows.append(parts[:7] if len(parts) >= 7 else parts + [""])
    df = pd.DataFrame(rows, columns=cols)
    df["dec_year"] = df["dec_year"].astype(float)
    df["ssn"] = df["ssn"].astype(float)
    df["date"] = pd.to_datetime(dict(year=df["year"].astype(int), month=df["month"].astype(int), day=15))
    return df[["date", "ssn"]].sort_values("date").reset_index(drop=True)


def ssn_for_dates(dates: pd.Series, ssn_monthly: pd.DataFrame) -> np.ndarray:
    # linear interpolation of monthly SSN onto arbitrary dates.
    # IMPORTANT: force both sides to the same datetime64 resolution before
    # casting to int64 -- pandas 2.x can hand back datetime64[us] from
    # pd.to_datetime(dict(...)) while a Timestamp+Timedelta sum stays [ns];
    # mixing the two silently misaligns the int64 epoch units (1000x off),
    # which makes np.interp saturate to a single boundary value for every
    # query (caught here: every sol was getting the same SSN, 2026-08's value).
    x = ssn_monthly["date"].values.astype("datetime64[ns]").astype("int64")
    y = ssn_monthly["ssn"].values
    xq = dates.values.astype("datetime64[ns]").astype("int64")
    return np.interp(xq, x, y)


def main():
    dose = load_dose()
    pressure = load_pressure()
    ssn_monthly = load_ssn_monthly()
    print(f"SSN monthly coverage: {ssn_monthly['date'].min().date()} .. {ssn_monthly['date'].max().date()}")

    max_sol = int(dose["sol"].max())
    full = pd.DataFrame({"sol": np.arange(0, max_sol + 1)})
    full = full.merge(dose, on="sol", how="left").merge(pressure, on="sol", how="left")
    full["date"] = LANDING_DATE + pd.to_timedelta(full["sol"] * SOL_TO_DAYS, unit="D")
    full["ssn"] = ssn_for_dates(full["date"], ssn_monthly)

    out_path = HERE / "mslrad_master_sol_series.csv"
    full.to_csv(out_path, index=False)
    print(f"saved {out_path} ({len(full)} rows, sol 0..{max_sol})")

    n_with_dose = full["dose_B_mean"].notna().sum() + full["dose_E_mean"].notna().sum()
    n_with_pressure = full["pressure_pa_mean"].notna().sum()
    print(f"sols with dose data: B={full['dose_B_mean'].notna().sum()}, E={full['dose_E_mean'].notna().sum()} "
          f"(of {max_sol+1})")
    print(f"sols with pressure data: {n_with_pressure} (of {max_sol+1})")
    print(full[["sol", "date", "dose_B_mean", "dose_E_mean", "pressure_pa_mean", "ssn"]].head())
    print(full[["sol", "date", "dose_B_mean", "dose_E_mean", "pressure_pa_mean", "ssn"]].tail())


if __name__ == "__main__":
    main()
