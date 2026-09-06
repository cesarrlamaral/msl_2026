"""
rad_extract_lib.py

Extracao de DOSIMETRY_TOTAL_DOSE_B/E dos RDR do MSL/RAD (PDS-PPI) via
HTTP Range requests direcionados por offset de byte, lidos do .LBL (detached
PDS3 label) de cada sol -- SEM baixar o .TXT inteiro (~65 MB/sol).

Cada sol tem um par (.LBL, .TXT) em:
  https://pds-ppi.igpp.ucla.edu/data/MSL-M-RAD-3-RDR-V1.0/DATA/SOL_XXXXX_YYYYY/

O .LBL contem, para cada observacao N (000, 002, 004, ...), um ponteiro:
  ^OBSnnn_TOT_DOSE_B_ELEMENT = ("ARQUIVO.TXT", OFFSET <BYTES>)
  OBJECT = OBSnnn_TOT_DOSE_B_ELEMENT
    START_BYTE = 1
    BYTES      = 8
  END_OBJECT

Posicao absoluta (1-indexed) do valor no .TXT = OFFSET + START_BYTE - 1.
Convertido para Range HTTP 0-indexed: start = OFFSET + START_BYTE - 2,
end = start + BYTES - 1.

Validado manualmente em 2026-09-05 contra SOL_02483_02579/RAD_RDR_2019_213_02_09_2483_V00:
  range 18891-18898 -> "11.39252" (TOT_DOSE_B, OBS000)
  range 18935-18942 -> "12.54804" (TOT_DOSE_E, OBS000)
Multi-range em uma unica requisicao HTTP confirmado funcional (resposta
multipart/byteranges, nginx no servidor pds-ppi.igpp.ucla.edu).
"""

import re
import time
import logging
from dataclasses import dataclass

import requests

BASE = "https://pds-ppi.igpp.ucla.edu/data/MSL-M-RAD-3-RDR-V1.0/DATA"
UA = "Mozilla/5.0 (research script; MSL RAD dose reanalysis; contact: cesarrlamaral@gmail.com)"

PTR_RE = re.compile(
    r'\^(OBS(\d+)_TOT_DOSE_([BE])_ELEMENT)\s*=\s*\n?\s*\("([^"]+)"\s*,\s*(\d+)\s*<BYTES>\)',
)
# after a pointer, the OBJECT block gives START_BYTE / BYTES for that same name
OBJ_RE_TMPL = r'OBJECT\s*=\s*{name}\s*\n(?:.*\n)*?\s*START_BYTE\s*=\s*(\d+)\s*\n\s*BYTES\s*=\s*(\d+)'

log = logging.getLogger("rad_extract")


@dataclass
class DoseObs:
    obs_index: int
    detector: str  # "B" or "E"
    range_start: int
    range_end: int
    value: float | None = None


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": UA})
    return s


def list_sol_range_folders(session: requests.Session) -> list[str]:
    """Return folder names like 'SOL_02483_02579' from the DATA/ listing."""
    r = session.get(BASE + "/", timeout=60)
    r.raise_for_status()
    names = sorted(set(re.findall(r'SOL_(\d{5}_\d{5})', r.text)))
    return [f"SOL_{n}" for n in names]


def list_sol_files(session: requests.Session, folder: str) -> list[str]:
    """Return sorted base names (without extension) of RAD_RDR_*_V*.TXT files in a folder."""
    r = session.get(f"{BASE}/{folder}/", timeout=60)
    r.raise_for_status()
    names = sorted(set(re.findall(r'(RAD_RDR_[A-Za-z0-9_]+_V\d+)\.TXT', r.text)))
    return names


def sol_number_from_name(basename: str) -> int:
    # RAD_RDR_2019_213_02_09_2483_V00 -> sol = 2483
    m = re.match(r'RAD_RDR_(?:C?\d+)_(\d+)_(\d+)_(\d+)_(\d+)_V\d+', basename)
    if not m:
        raise ValueError(f"unexpected filename pattern: {basename}")
    return int(m.group(4))


def fetch_dose_offsets(session: requests.Session, folder: str, basename: str) -> list[DoseObs]:
    """Download the .LBL for one sol and parse byte ranges for every TOT_DOSE_B/E element."""
    url = f"{BASE}/{folder}/{basename}.LBL"
    r = session.get(url, timeout=60)
    r.raise_for_status()
    text = r.text

    obs_list: list[DoseObs] = []
    for m in PTR_RE.finditer(text):
        obj_name, obs_idx, det, fname, offset = m.groups()
        offset = int(offset)
        obs_idx = int(obs_idx)
        obj_re = re.compile(OBJ_RE_TMPL.format(name=re.escape(obj_name)))
        om = obj_re.search(text, m.end())
        if not om:
            log.warning("no OBJECT block found for %s in %s/%s", obj_name, folder, basename)
            continue
        start_byte, n_bytes = int(om.group(1)), int(om.group(2))
        abs0 = offset + start_byte - 2  # 0-indexed HTTP range start
        obs_list.append(DoseObs(obs_idx, det, abs0, abs0 + n_bytes - 1))
    return obs_list


def fetch_dose_values(session: requests.Session, folder: str, basename: str,
                       obs_list: list[DoseObs], chunk: int = 120) -> None:
    """Fill in .value for each DoseObs via one (or a few, chunked) multi-range GET(s)."""
    url = f"{BASE}/{folder}/{basename}.TXT"
    by_range = {(o.range_start, o.range_end): o for o in obs_list}
    ranges = list(by_range.keys())

    for i in range(0, len(ranges), chunk):
        part = ranges[i:i + chunk]
        range_header = "bytes=" + ",".join(f"{a}-{b}" for a, b in part)
        r = session.get(url, headers={"Range": range_header}, timeout=90)
        r.raise_for_status()
        ctype = r.headers.get("Content-Type", "")
        if "multipart/byteranges" not in ctype:
            # server collapsed to a single part (e.g. chunk of 1) -- handle directly
            if len(part) == 1:
                by_range[part[0]].value = float(r.content.decode().strip())
                continue
            raise RuntimeError(f"expected multipart response, got {ctype}")
        boundary = ctype.split("boundary=")[1].strip()
        _parse_multipart(r.content, boundary, by_range)


def _parse_multipart(body: bytes, boundary: str, by_range: dict) -> None:
    marker = ("--" + boundary).encode()
    parts = body.split(marker)
    cr_re = re.compile(rb"Content-range:\s*bytes\s*(\d+)-(\d+)/\d+", re.IGNORECASE)
    for part in parts:
        m = cr_re.search(part)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2))
        payload = part.split(b"\r\n\r\n", 1)[-1].strip(b"\r\n-")
        key = (a, b)
        if key in by_range:
            try:
                by_range[key].value = float(payload.decode().strip())
            except ValueError:
                log.warning("could not parse value %r for range %s", payload, key)


def process_folder_to_csv(session: requests.Session, folder: str, out_path, n_workers: int = 6) -> None:
    """Extract every sol file in a folder into out_path (CSV: sol,obs_index,detector,value).

    IMPORTANT: some sol-range folders have MULTIPLE RDR files for the same sol
    (different DSN downlink passes covering the same sol; seen e.g. in
    SOL_01649_01772, up to 17 files for one sol). obs_index alone is only
    unique WITHIN a single file, so it cannot be used as-is to key rows across
    files of the same sol -- doing so silently collapses distinct
    observations. Here, files sharing a sol are sorted chronologically (their
    names encode year/day/hour/minute) and each gets a fixed slot: global
    obs_index = file_slot*1000 + in_file_obs_index.

    Resumability is tracked per SOURCE FILE (a sidecar .manifest.txt next to
    the CSV, one processed basename per line) rather than per sol, precisely
    because a sol can need more than one file before it is complete.
    """
    import csv as _csv
    from pathlib import Path as _Path
    from collections import defaultdict as _defaultdict
    from concurrent.futures import ThreadPoolExecutor as _TPE, as_completed as _as_completed

    out_path = _Path(out_path)
    manifest_path = out_path.with_suffix(".manifest.txt")

    files = list_sol_files(session, folder)
    by_sol: dict[int, list[str]] = _defaultdict(list)
    for b in files:
        by_sol[sol_number_from_name(b)].append(b)
    file_slot = {b: i for blist in by_sol.values() for i, b in enumerate(sorted(blist))}

    done_files = set(manifest_path.read_text().splitlines()) if manifest_path.exists() else set()
    todo = [b for b in files if b not in done_files]
    log.info("[%s] %d files total, %d already done, %d to fetch", folder, len(files), len(done_files), len(todo))
    if not todo:
        return

    write_header = not out_path.exists()
    with open(out_path, "a", newline="") as f, open(manifest_path, "a") as mf:
        writer = _csv.writer(f)
        if write_header:
            writer.writerow(["sol", "obs_index", "detector", "value"])

        with _TPE(max_workers=n_workers) as pool:
            futures = {pool.submit(_safe_extract_sol, session, folder, b): b for b in todo}
            n_ok, n_fail = 0, 0
            for fut in _as_completed(futures):
                basename = futures[fut]
                result = fut.result()
                if result is None:
                    n_fail += 1
                    continue
                sol, obs = result
                slot = file_slot[basename]
                for o in sorted(obs, key=lambda o: (o.obs_index, o.detector)):
                    global_idx = slot * 1000 + o.obs_index
                    writer.writerow([sol, global_idx, o.detector, o.value])
                f.flush()
                mf.write(basename + "\n")
                mf.flush()
                n_ok += 1
                if n_ok % 25 == 0:
                    log.info("[%s] progress: %d/%d ok, %d failed", folder, n_ok, len(todo), n_fail)
    log.info("[%s] DONE: %d ok, %d failed", folder, n_ok, n_fail)


def _safe_extract_sol(session, folder, basename):
    try:
        return extract_sol(session, folder, basename)
    except Exception as e:
        log.warning("[%s/%s] extraction failed: %s", folder, basename, e)
        return None


def mark_folder_fully_done(folder: str, out_dir, session: requests.Session) -> None:
    """For a folder whose CSV is already known-complete and duplicate-free
    (verified separately), write a manifest listing all its basenames so
    process_folder_to_csv treats it as done and does not reprocess/duplicate it."""
    from pathlib import Path as _Path
    out_path = _Path(out_dir) / f"{folder}.csv"
    manifest_path = out_path.with_suffix(".manifest.txt")
    if manifest_path.exists():
        return
    files = list_sol_files(session, folder)
    manifest_path.write_text("\n".join(files) + "\n")


def extract_sol(session: requests.Session, folder: str, basename: str,
                 retries: int = 3, backoff: float = 2.0) -> tuple[int, list[DoseObs]]:
    sol = sol_number_from_name(basename)
    last_exc = None
    for attempt in range(retries):
        try:
            obs = fetch_dose_offsets(session, folder, basename)
            if obs:
                fetch_dose_values(session, folder, basename, obs)
            return sol, obs
        except (requests.RequestException, RuntimeError) as e:
            last_exc = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"failed to extract {folder}/{basename} after {retries} tries: {last_exc}")
