"""
04_validate_extraction.py

Valida o pipeline de extracao original (mslrad_v1.py, heuristica de "2 linhas
apos o header de texto") contra o metodo novo, byte-exato, guiado pelos
ponteiros START_BYTE/BYTES do .LBL do PDS-PPI (rad_extract_lib.py).

Estrategia: re-extrai, sol a sol, os primeiros N sols da janela SOL_00000_00089
(a mais antiga que ja temos localmente processada) e compara a sequencia de
valores contra o inicio de results_detector_B.txt / results_detector_E.txt
daquela janela (MSL/2026/MSLRADv1 RESULTS/SOL 00000_00089/).
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from rad_extract_lib import make_session, list_sol_files, extract_sol  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("validate")

FOLDER = "SOL_00000_00089"
LOCAL_DIR = Path(r"C:\Users\Usuário\Desktop\Projetos_Claude\MSL\2026\MSLRADv1 RESULTS\SOL 00000_00089")
N_SOLS_TO_CHECK = 999  # 999 = janela inteira (78 arquivos)


def load_local_series(path: Path) -> list[float]:
    with open(path, "r") as f:
        return [float(line.strip()) for line in f if line.strip()]


def main():
    session = make_session()
    all_files = list_sol_files(session, FOLDER)
    log.info("folder %s has %d sol files total; checking first %d", FOLDER, len(all_files), N_SOLS_TO_CHECK)

    subset = all_files[:N_SOLS_TO_CHECK]
    new_B, new_E = [], []
    per_sol_counts = []
    for basename in subset:
        sol, obs = extract_sol(session, FOLDER, basename)
        b_vals = [o.value for o in sorted(obs, key=lambda o: o.obs_index) if o.detector == "B"]
        e_vals = [o.value for o in sorted(obs, key=lambda o: o.obs_index) if o.detector == "E"]
        new_B.extend(b_vals)
        new_E.extend(e_vals)
        per_sol_counts.append((sol, len(b_vals), len(e_vals)))
        log.info("sol %d: %d B obs, %d E obs", sol, len(b_vals), len(e_vals))

    local_B = load_local_series(LOCAL_DIR / "results_detector_B.txt")
    local_E = load_local_series(LOCAL_DIR / "results_detector_E.txt")

    log.info("=== COMPARISON (detector B) ===")
    log.info("new method produced %d values; comparing against first %d of local (%d total)",
              len(new_B), len(new_B), len(local_B))
    n = min(len(new_B), len(local_B))
    mismatches_B = [(i, new_B[i], local_B[i]) for i in range(n) if abs(new_B[i] - local_B[i]) > 1e-4]
    log.info("B: %d/%d compared values match within 1e-4; %d mismatches", n - len(mismatches_B), n, len(mismatches_B))
    for i, nv, lv in mismatches_B[:10]:
        log.info("  mismatch idx %d: new=%s local=%s", i, nv, lv)

    log.info("=== COMPARISON (detector E) ===")
    n2 = min(len(new_E), len(local_E))
    mismatches_E = [(i, new_E[i], local_E[i]) for i in range(n2) if abs(new_E[i] - local_E[i]) > 1e-4]
    log.info("E: %d/%d compared values match within 1e-4; %d mismatches", n2 - len(mismatches_E), n2, len(mismatches_E))
    for i, nv, lv in mismatches_E[:10]:
        log.info("  mismatch idx %d: new=%s local=%s", i, nv, lv)

    log.info("per-sol observation counts: %s", per_sol_counts)
    log.info("DONE. new_B total=%d new_E total=%d local_B total=%d local_E total=%d",
              len(new_B), len(new_E), len(local_B), len(local_E))


if __name__ == "__main__":
    main()
