"""
06_extract_rems_pressure.py

Extrai pressao do REMS (produto RMD) para TODA a faixa de sols coberta pelo
paper + a extensao (SOL 1 ate a mais recente disponivel), amostrando
n_samples pontos por sol via HTTP Range request (ver rems_extract_lib.py).

Isso cobre tanto os sols originais do paper (1-2482, para permitir refazer
o teste diurno dose-vs-pressao com dado real em vez de so citar Harri 2014)
quanto a extensao (2483-4843, para o teste de regressao dose-vs-pressao no
residuo apos remover a tendencia solar -- ponto (1) do plano).

Saida: um CSV por janela SOL_XXXXX_YYYYY em pressure_series/, colunas:
sol,row_index,timestamp,ltst,pressure_pa. Resumivel.
"""

import csv
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from rems_extract_lib import make_session, list_sol_range_folders, list_sols_in_folder, extract_sol_pressure  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("rems_extend")

OUT_DIR = Path(__file__).parent / "pressure_series"
OUT_DIR.mkdir(exist_ok=True)
N_SAMPLES_PER_SOL = 48
N_WORKERS = 8


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
    sols = list_sols_in_folder(session, folder)
    done_sols = load_done_sols(out_path)
    todo = [s for s in sols if s not in done_sols]
    log.info("[%s] %d sols total, %d already done, %d to fetch", folder, len(sols), len(done_sols), len(todo))

    write_header = not out_path.exists()
    n_ok, n_missing, n_fail = 0, 0, 0

    def _safe(sol):
        try:
            return sol, extract_sol_pressure(session, folder, sol, n_samples=N_SAMPLES_PER_SOL), None
        except Exception as e:
            return sol, None, e

    with open(out_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["sol", "row_index", "timestamp", "ltst", "pressure_pa"])
        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = [pool.submit(_safe, sol) for sol in todo]
            done_count = 0
            for fut in as_completed(futures):
                sol, samples, err = fut.result()
                done_count += 1
                if err is not None:
                    log.warning("[%s] sol %d failed: %s", folder, sol, err)
                    n_fail += 1
                    continue
                if samples is None:
                    n_missing += 1
                    continue
                for s in samples:
                    writer.writerow([sol, s.row_index, s.timestamp, s.ltst, s.pressure])
                f.flush()
                n_ok += 1
                if done_count % 25 == 0:
                    log.info("[%s] progress: %d/%d done (%d ok, %d missing, %d failed)",
                              folder, done_count, len(todo), n_ok, n_missing, n_fail)
    log.info("[%s] DONE: %d ok, %d missing, %d failed", folder, n_ok, n_missing, n_fail)


def main():
    session = make_session()
    all_folders = list_sol_range_folders(session)
    log.info("REMS archive has %d sol-range folders total: %s ... %s", len(all_folders), all_folders[0], all_folders[-1])

    t0 = time.time()
    for folder in all_folders:
        process_folder(session, folder)
    log.info("ALL DONE in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
