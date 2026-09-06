"""
rems_extract_lib.py

Extracao da pressao do REMS (MSL) via HTTP Range requests, para regredir a
dose do RAD contra pressao real (em vez de so citar Harri et al. 2014 em
prosa) -- ponto (1) do plano de fortalecimento do paper.

Fonte: https://atmos.nmsu.edu/PDS/data/mslrem_1001/DATA/SOL_XXXXX_YYYYY/SOLnnnnn/
Produto usado: RMD ("MOD RDR" -- pressao unica corrigida/modelada, Pascal),
NAO o RNV (que so tem os dois sensores BAROCAP brutos separados).

Schema PDS4 confirmado constante (amostrado em SOL 2483 e SOL 4710):
  record_length = 372 bytes, offset = 0 (tabela comeca no byte 0 do .TAB)
  TIMESTAMP:  byte 1,   len 10  (segundos desde 1-Jan-2000 12:00)
  LTST:       byte 34,  len 14  (Local True Solar Time)
  PRESSURE:   byte 345, len 7   (Pascal)
Byte offsets no .xml sao 1-indexed; para HTTP Range (0-indexed) subtrai-se 1.

Como o schema e fixo, NAO e preciso baixar o .xml de cada sol -- so a
listagem HTML do diretorio do sol (pra descobrir o nome do .TAB, que tem um
prefixo de clock de espaconave imprevisivel) e o tamanho do arquivo (pra
saber quantos records tem).
"""

import re
import time
import logging
from dataclasses import dataclass

import requests

BASE = "https://atmos.nmsu.edu/PDS/data/mslrem_1001/DATA"
UA = "Mozilla/5.0 (research script; MSL REMS pressure reanalysis; contact: cesarrlamaral@gmail.com)"

RECORD_LENGTH = 372
FIELDS = {
    "TIMESTAMP": (1, 10),
    "LTST": (34, 14),
    "PRESSURE": (345, 7),
}

log = logging.getLogger("rems_extract")


@dataclass
class PressureSample:
    row_index: int
    timestamp: str | None = None
    ltst: str | None = None
    pressure: float | None = None


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": UA})
    return s


def list_sol_range_folders(session: requests.Session) -> list[str]:
    r = session.get(BASE + "/", timeout=60)
    r.raise_for_status()
    names = sorted(set(re.findall(r'SOL_(\d{5}_\d{5})', r.text)))
    return [f"SOL_{n}" for n in names]


def list_sols_in_folder(session: requests.Session, folder: str) -> list[int]:
    r = session.get(f"{BASE}/{folder}/", timeout=60)
    r.raise_for_status()
    return sorted(int(s) for s in set(re.findall(r'href="SOL(\d+)/"', r.text)))


def find_rmd_file(session: requests.Session, folder: str, sol: int) -> tuple[str, int] | None:
    """Return (tab_basename_without_ext, file_size_bytes) for the RMD product of a sol, or None."""
    sol_dir = f"SOL{sol:05d}"
    r = session.get(f"{BASE}/{folder}/{sol_dir}/", timeout=60)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    m = re.search(r'href="([^"]*RMD[^"]*)\.TAB"', r.text)
    if not m:
        return None
    basename = m.group(1)
    size_m = re.search(re.escape(f'{basename}.TAB') + r'.*?(\d+)</td>\s*<td[^>]*>\s*(\d+)</td>', r.text)
    return basename, sol_dir


def get_file_size(session: requests.Session, folder: str, sol_dir: str, basename: str) -> int:
    url = f"{BASE}/{folder}/{sol_dir}/{basename}.TAB"
    r = session.head(url, timeout=60)
    r.raise_for_status()
    return int(r.headers["Content-Length"])


def build_sample_rows(n_records: int, n_samples: int = 48) -> list[int]:
    if n_records <= n_samples:
        return list(range(n_records))
    step = n_records / n_samples
    return sorted(set(int(i * step) for i in range(n_samples)))


def row_byte_ranges(row: int) -> dict[str, tuple[int, int]]:
    ranges = {}
    base = row * RECORD_LENGTH
    for field, (byte1, length) in FIELDS.items():
        start0 = base + byte1 - 1  # 1-indexed -> 0-indexed
        ranges[field] = (start0, start0 + length - 1)
    return ranges


def fetch_samples(session: requests.Session, folder: str, sol_dir: str, basename: str,
                   rows: list[int], chunk: int = 100) -> list[PressureSample]:
    url = f"{BASE}/{folder}/{sol_dir}/{basename}.TAB"
    samples = {row: PressureSample(row_index=row) for row in rows}

    flat_ranges = []
    for row in rows:
        for field, rng in row_byte_ranges(row).items():
            flat_ranges.append((rng, row, field))

    for i in range(0, len(flat_ranges), chunk):
        part = flat_ranges[i:i + chunk]
        range_header = "bytes=" + ",".join(f"{a}-{b}" for (a, b), _, _ in part)
        r = session.get(url, headers={"Range": range_header}, timeout=90)
        r.raise_for_status()
        ctype = r.headers.get("Content-Type", "")
        if "multipart/byteranges" not in ctype:
            if len(part) == 1:
                (a, b), row, field = part[0]
                _assign(samples[row], field, r.content)
                continue
            raise RuntimeError(f"expected multipart, got {ctype}")
        boundary = ctype.split("boundary=")[1].strip()
        _parse_multipart(r.content, boundary, part, samples)

    return [samples[row] for row in rows]


def _assign(sample: PressureSample, field: str, raw: bytes) -> None:
    text = raw.decode(errors="replace").strip()
    if field == "PRESSURE":
        try:
            sample.pressure = float(text)
        except ValueError:
            sample.pressure = None
    elif field == "TIMESTAMP":
        sample.timestamp = text
    elif field == "LTST":
        sample.ltst = text


def _parse_multipart(body: bytes, boundary: str, part_defs, samples: dict) -> None:
    key_to_target = {(a, b): (row, field) for (a, b), row, field in part_defs}
    marker = ("--" + boundary).encode()
    cr_re = re.compile(rb"Content-range:\s*bytes\s*(\d+)-(\d+)/\d+", re.IGNORECASE)
    for chunk in body.split(marker):
        m = cr_re.search(chunk)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2))
        payload = chunk.split(b"\r\n\r\n", 1)[-1].strip(b"\r\n-")
        key = (a, b)
        if key in key_to_target:
            row, field = key_to_target[key]
            _assign(samples[row], field, payload)


def extract_sol_pressure(session: requests.Session, folder: str, sol: int,
                          n_samples: int = 48, retries: int = 3) -> list[PressureSample] | None:
    sol_dir = f"SOL{sol:05d}"
    last_exc = None
    for attempt in range(retries):
        try:
            r = session.get(f"{BASE}/{folder}/{sol_dir}/", timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            m = re.search(r'href="([^"]*RMD[^"]*)\.TAB"', r.text)
            if not m:
                return None
            basename = m.group(1)
            size = get_file_size(session, folder, sol_dir, basename)
            n_records = size // RECORD_LENGTH
            rows = build_sample_rows(n_records, n_samples)
            return fetch_samples(session, folder, sol_dir, basename, rows)
        except (requests.RequestException, RuntimeError) as e:
            last_exc = e
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"failed to extract pressure for {folder}/sol {sol}: {last_exc}")
