"""
07_backfill_original_history.py

Reextrai os sols 1-2482 (todas as janelas ja usadas no paper original) com o
mesmo metodo byte-exato usado em 05_extend_series.py, para dar identidade de
sol/observacao a TODA a serie -- nao so a parte estendida.

Motivacao: o pipeline de 2020 (mslrad_v1.py) so concatenava valores, sem
guardar a que sol cada observacao pertencia. Isso forcava a decomposicao/
wavelet original a tratar "indice da observacao" como proxy de tempo,
assumindo cadencia uniforme (que na pratica varia de 13 a 88 obs/sol -- ver
04_validate_extraction.py). Com sol real por observacao, sol 1 a 4843 vira
uma serie temporal unica e propriamente indexada.

Ja validado byte-a-byte contra os results_detector_B/E.txt originais em
04_validate_extraction.py (3763/3763 valores identicos).

BUG corrigido nesta versao (2026-09-06): algumas janelas tem MULTIPLOS
arquivos de sol para o MESMO sol (passes de downlink diferentes -- ate 17
arquivos para um unico sol em SOL_01649_01772). A primeira versao deste
script usava so "sol" como chave de resumo/deduplicacao, o que descartava
observacoes reais de arquivos extras do mesmo sol como se fossem duplicatas.
Corrigido em rad_extract_lib.process_folder_to_csv (ver docstring la):
obs_index agora e unico por (sol, slot_do_arquivo), e o resumo e rastreado
por ARQUIVO fonte (manifest), nao por sol. As 15 janelas ja processadas
corretamente antes da correcao (confirmado 0 sols duplicados nelas) sao
marcadas como concluidas sem reprocessar; as 4 janelas afetadas
(SOL_00180_00269, SOL_00939_01062, SOL_01294_01417, SOL_01649_01772) tiveram
o CSV anterior apagado e sao reextraidas do zero.
"""

import logging
import time
from pathlib import Path

from rad_extract_lib import make_session, list_sol_range_folders, process_folder_to_csv, mark_folder_fully_done

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("backfill")

OUT_DIR = Path(__file__).parent / "extended_series"
OUT_DIR.mkdir(exist_ok=True)

LAST_ORIGINAL_FOLDER = "SOL_02359_02482"

# Janelas confirmadas (diagnostico de 2026-09-06, ver EDITORIAL_NOTES.md) como
# ja extraidas ANTES da correcao e SEM nenhum sol com >1 arquivo -- portanto
# corretas como estao. Marcadas como concluidas via manifest, nao reprocessadas.
CLEAN_ALREADY_DONE = [
    "SOL_00000_00089", "SOL_00090_00179", "SOL_00270_00359", "SOL_00360_00449",
    "SOL_00450_00583", "SOL_00584_00707", "SOL_00708_00804", "SOL_00805_00938",
    "SOL_01063_01159", "SOL_01160_01293", "SOL_01418_01514", "SOL_01515_01648",
]


def main():
    session = make_session()
    all_folders = list_sol_range_folders(session)
    end_idx = all_folders.index(LAST_ORIGINAL_FOLDER) + 1
    original_folders = all_folders[:end_idx]
    log.info("original range: %d folders: %s ... %s", len(original_folders), original_folders[0], original_folders[-1])

    for folder in CLEAN_ALREADY_DONE:
        mark_folder_fully_done(folder, OUT_DIR, session)
        log.info("[%s] marked as already complete (verified duplicate-free)", folder)

    t0 = time.time()
    for folder in original_folders:
        process_folder_to_csv(session, folder, OUT_DIR / f"{folder}.csv")
    log.info("ALL DONE in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
