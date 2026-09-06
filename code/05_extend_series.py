"""
05_extend_series.py

Estende a serie de dose do RAD (detectores B e E) alem do SOL 2482 usando o
mesmo metodo validado em 04_validate_extraction.py (extracao byte-exata via
.LBL do PDS-PPI, sem baixar os .TXT de 65 MB/sol inteiros).

Cobre todas as janelas SOL_XXXXX_YYYYY publicadas apos a ultima usada no
paper (SOL_02359_02482), ate a mais recente disponivel no archive.

Saida: um CSV por janela, com colunas sol,obs_index,detector,value -- preserva
a identidade de sol e observacao (o pipeline original de 2020 nao guardava
isso, so concatenava valores). Resumivel: pula janelas/sols ja completos.
"""

import csv
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from rad_extract_lib import make_session, list_sol_range_folders, list_sol_files, extract_sol  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("extend")

OUT_DIR = Path(__file__).parent / "extended_series"
OUT_DIR.mkdir(exist_ok=True)

LAST_ORIGINAL_FOLDER = "SOL_02359_02482"  # ultima janela ja usada no paper
N_WORKERS = 6


def folder_csv_path(folder: str) -> Path:
    return OUT_DIR / f"{folder}.csv"


def load_done_sols(csv_path: Path) -> set[int]:
    if not csv_path.exists():
        return set()
    done = set()
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            done.add(int(row["sol"]))
    return done


def process_folder(session, folder: str) -> None:
    out_path = folder_csv_path(folder)
    files = list_sol_files(session, folder)
    done_sols = load_done_sols(out_path)
    todo = [b for b in files if _sol_of(b) not in done_sols]
    log.info("[%s] %d files total, %d already done, %d to fetch", folder, len(files), len(done_sols), len(todo))

    write_header = not out_path.exists()
    with open(out_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["sol", "obs_index", "detector", "value"])

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(_safe_extract, session, folder, b): b for b in todo}
            n_ok, n_fail = 0, 0
            for fut in as_completed(futures):
                basename = futures[fut]
                result = fut.result()
                if result is None:
                    n_fail += 1
                    continue
                sol, obs = result
                for o in sorted(obs, key=lambda o: (o.obs_index, o.detector)):
                    writer.writerow([sol, o.obs_index, o.detector, o.value])
                f.flush()
                n_ok += 1
                if n_ok % 25 == 0:
                    log.info("[%s] progress: %d/%d ok, %d failed", folder, n_ok, len(todo), n_fail)
    log.info("[%s] DONE: %d ok, %d failed", folder, n_ok, n_fail)


def _sol_of(basename: str) -> int:
    from rad_extract_lib import sol_number_from_name
    return sol_number_from_name(basename)


def _safe_extract(session, folder, basename):
    try:
        return extract_sol(session, folder, basename)
    except Exception as e:
        log.warning("[%s/%s] extraction failed: %s", folder, basename, e)
        return None


def main():
    session = make_session()
    all_folders = list_sol_range_folders(session)
    log.info("archive has %d sol-range folders total", len(all_folders))

    start_idx = all_folders.index(LAST_ORIGINAL_FOLDER) + 1
    new_folders = all_folders[start_idx:]
    log.info("extending with %d new folders: %s ... %s", len(new_folders), new_folders[0], new_folders[-1])

    t0 = time.time()
    for folder in new_folders:
        process_folder(session, folder)
    log.info("ALL DONE in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
