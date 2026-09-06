# msl_2026

Extraction pipeline, reanalysis code, and manuscript for a paper analyzing the
Mars Science Laboratory (MSL) Curiosity rover's Radiation Assessment Detector
(RAD) dose-rate record, SOL 1-4843 (Aug 2012 - Mar 2026).

## Contents

- `manuscript/` — the paper (`paper_msl_2026.docx`), editorial notes
  documenting the analysis process (`EDITORIAL_NOTES.md`), and all figures
  embedded in the manuscript (`figuras/`).
- `code/` — the Python extraction and analysis pipeline, numbered in the
  order it was run:
  - `rad_extract_lib.py`, `rems_extract_lib.py` — byte-offset HTTP Range
    extraction from the PDS archives (RAD: PDS-PPI node; REMS: PDS
    Atmospheres node), without downloading full per-sol RDR files.
  - `01`-`04` — original-record decomposition/wavelet/t-test reanalysis
    (R-to-Python migration) and byte-exact validation of the extraction
    against the PDS label pointers.
  - `05`-`08` — extending the record to SOL 4843 and consolidating dose,
    pressure, and solar-activity data into one per-sol series.
  - `09`-`10` — extended decomposition/wavelet analysis and the continuous
    solar-activity/pressure regression (Section 3.3-3.4 of the paper).
  - `11`-`12` — physical GCR proxy (Oulu neutron monitor), lag/hysteresis
    test, and Lomb-Scargle periodogram (Section 3.5-3.6).
  - `data/` — the key derived datasets these scripts produce and consume:
    the consolidated per-sol series (`mslrad_master_sol_series.csv`), the
    Oulu neutron monitor daily count rate (`nm_data/oulu_daily.csv`), and
    small result-summary CSVs.

## Data sources

- RAD and REMS raw data: NASA Planetary Data System (RAD: Planetary Plasma
  Interactions node, https://pds-ppi.igpp.ucla.edu; REMS: Atmospheres node,
  https://atmos.nmsu.edu/PDS/data/mslrem_1001/).
- Neutron monitor count rate: Oulu station, via the Neutron Monitor Database
  (NMDB, https://www.nmdb.eu).
- Observed sunspot number: SILSO World Data Center, Royal Observatory of
  Belgium.

Full citations for all data sources and methods are in the manuscript's
References and Data and Code Availability sections.

## Not included

Bibliography PDFs (copyrighted journal articles cited in the manuscript) and
large intermediate per-observation exports are not included in this
repository; the consolidated per-sol dataset in `code/data/` is sufficient to
reproduce the paper's figures and tables.
