#!/usr/bin/env python3
"""
Config-driven mitochondrial variant liftover + codon-match + tRNA-match pipeline.

Design:
  1. Build pairwise alignment and species->human position map.
  2. Liftover-only VCF/COV. Raw lifted VCF keeps original species coordinates in INFO.
  3. Optionally annotate codon-match. Annotation-only, no filtering.
  4. Optionally build tRNAscan-SE position indexes and annotate tRNA-match. Annotation-only.
  5. Optionally apply a final INFO-based filter.

Only Python standard library is required. External tools are needed only for optional steps:
  - mafft/muscle/prank for alignment
  - tRNAscan-SE for tRNA discovery
"""
from __future__ import annotations

import argparse
import bisect
import configparser
import csv
import gzip
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

VALID = set("ACGT")
NA = "NA"

# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------

def die(msg: str, code: int = 2) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)


def log(msg: str) -> None:
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"[WARN] {msg}", file=sys.stderr)


def safe_div(num, den):
    if den is None or den == 0:
        return "NA"
    return f"{num / den:.6f}"


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y", "on"}


def sanitize_info_value(x) -> str:
    s = str(x)
    s = s.replace(";", ",").replace("\t", "_").replace(" ", "_")
    if s == "":
        return NA
    return s


def parse_optional_int(x) -> Optional[int]:
    if x is None:
        return None
    s = str(x).strip()
    if s in {"", ".", NA}:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def parse_optional_float(x) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s in {"", ".", NA}:
        return None
    try:
        return float(s)
    except Exception:
        return None


def strip_vcf_suffix(path: str | Path) -> str:
    name = Path(path).name
    for suf in [".vcf.gz", ".vcf", ".bcf", ".gz"]:
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def write_tsv(path: str | Path, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    with open_text(path, "wt") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, NA) for k in fieldnames})


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

class PipeConfig:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        if not self.path.exists():
            die(f"Config not found: {self.path}")
        self.cp = configparser.ConfigParser(interpolation=None)
        self.cp.optionxform = str.lower
        self.cp.read(self.path)

    def has(self, section: str, key: str) -> bool:
        return self.cp.has_section(section) and self.cp.has_option(section, key)

    def get(self, section: str, key: str, default: Optional[str] = None, **ctx) -> str:
        if self.has(section, key):
            val = self.cp.get(section, key).strip()
        elif default is not None:
            val = str(default)
        else:
            die(f"Missing config value [{section}] {key}")
        return self.expand(val, **ctx)

    def getint(self, section: str, key: str, default: Optional[int] = None, **ctx) -> int:
        return int(self.get(section, key, str(default) if default is not None else None, **ctx))

    def getbool(self, section: str, key: str, default: bool = False, **ctx) -> bool:
        return parse_bool(self.get(section, key, "1" if default else "0", **ctx))

    def expand(self, val: str, **ctx) -> str:
        # Build a context with all path/settings values. Apply formatting a few times to resolve {outdir}-style paths.
        base: Dict[str, str] = {}
        for sec in self.cp.sections():
            for k, v in self.cp.items(sec):
                base[k] = v.strip()
        base.update({k: str(v) for k, v in ctx.items() if v is not None})
        out = val
        for _ in range(8):
            try:
                new = out.format(**base)
            except KeyError:
                break
            if new == out:
                break
            out = new
            base.update({k: out if k == "_" else v for k, v in base.items()})
        return os.path.expandvars(os.path.expanduser(out))

    def path_get(self, key: str, default: Optional[str] = None, **ctx) -> str:
        return self.get("paths", key, default, **ctx)

    def setting(self, key: str, default: str = "", **ctx) -> str:
        return self.get("settings", key, default, **ctx)

    def setting_bool(self, key: str, default: bool = False, **ctx) -> bool:
        return self.getbool("settings", key, default, **ctx)

    def setting_int(self, key: str, default: int = 0, **ctx) -> int:
        return self.getint("settings", key, default, **ctx)


def make_outdirs(cfg: PipeConfig) -> Dict[str, Path]:
    outdir = ensure_dir(cfg.path_get("outdir"))
    d = {
        "outdir": outdir,
        "tmp": ensure_dir(cfg.path_get("tmp_dir", "{outdir}/tmp")),
        "reports": ensure_dir(cfg.path_get("reports_dir", "{outdir}/reports")),
        "maps": ensure_dir(cfg.path_get("maps_dir", "{outdir}/maps")),
        "alignments": ensure_dir(cfg.path_get("alignments_dir", "{outdir}/alignments")),
        "vcf_raw": ensure_dir(cfg.path_get("vcf_lifted_raw_dir", "{outdir}/vcf_lifted_raw")),
        "vcf_codon": ensure_dir(cfg.path_get("vcf_codon_dir", "{outdir}/vcf_annotated_codon")),
        "vcf_trna": ensure_dir(cfg.path_get("vcf_trna_dir", "{outdir}/vcf_annotated_trna")),
        "vcf_final": ensure_dir(cfg.path_get("vcf_final_dir", "{outdir}/vcf_final")),
        "cov": ensure_dir(cfg.path_get("cov_lifted_dir", "{outdir}/cov_lifted")),
        "trna_index": ensure_dir(cfg.path_get("trna_index_dir", "{outdir}/trna_index")),
        "debug": ensure_dir(cfg.path_get("debug_dir", "{outdir}/debug")),
    }
    return d


# -----------------------------------------------------------------------------
# FASTA and sample metadata
# -----------------------------------------------------------------------------

def read_fasta_all(path: str | Path) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, str]] = []
    name = None
    seq: List[str] = []
    with open_text(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq).upper()))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            recs.append((name, "".join(seq).upper()))
    return recs


def first_fasta_record(path: str | Path) -> Tuple[str, str]:
    recs = read_fasta_all(path)
    if not recs:
        die(f"No FASTA records in {path}")
    if len(recs) > 1:
        warn(f"{path} has {len(recs)} records; using the first one: {recs[0][0]}")
    return recs[0]


def load_sample_ref_map(file_path: str | Path) -> Dict[str, str]:
    # Compatible with the user's current file convention: col1 accession, col3 species/sample key.
    mp: Dict[str, str] = {}
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 3:
                mp[p[2]] = p[2]
    return mp


def load_species_to_accessions(file_path: str | Path) -> Dict[str, List[str]]:
    mp: Dict[str, List[str]] = defaultdict(list)
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 3:
                mp[p[2]].append(p[0])
    return {k: sorted(set(v)) for k, v in mp.items()}


def load_rotate_info(file_path: str | Path) -> Dict[str, Tuple[int, int]]:
    # Returns key -> (length, init_position). This matches the user's current script.
    info: Dict[str, Tuple[int, int]] = {}
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 3:
                info[p[0]] = (int(p[1]), int(p[2]))
    return info


def discover_samples(cfg: PipeConfig) -> List[str]:
    sample_list = cfg.path_get("sample_list_file", "")
    if sample_list and Path(sample_list).exists():
        out = []
        with open(sample_list) as f:
            for line in f:
                s = line.strip().split()
                if s and not s[0].startswith("#"):
                    out.append(s[0])
        return out
    primate_dir = cfg.path_get("primate_dir")
    hits = sorted(glob.glob(os.path.join(primate_dir, "*.fa")) + glob.glob(os.path.join(primate_dir, "*.fasta")))
    return [Path(h).stem for h in hits]


def sample_fasta(cfg: PipeConfig, sample: str) -> str:
    template = cfg.path_get("sample_fasta_template", "", sample=sample)
    if template and Path(template).exists():
        return template
    primate_dir = cfg.path_get("primate_dir")
    for ext in [".fa", ".fasta"]:
        p = Path(primate_dir) / f"{sample}{ext}"
        if p.exists():
            return str(p)
    die(f"Cannot find FASTA for sample={sample} in {primate_dir}")
    return ""


def sample_length_from_fai_or_fasta(cfg: PipeConfig, sample: str, fasta: str) -> int:
    fai_dir = cfg.path_get("fai_dir", "{primate_dir}", sample=sample)
    candidates = [
        Path(fai_dir) / f"{sample}.fa.fai",
        Path(fai_dir) / f"{sample}.fasta.fai",
        Path(str(fasta) + ".fai"),
    ]
    for fai in candidates:
        if fai.exists():
            with open(fai) as f:
                first = f.readline().rstrip("\n").split("\t")
                if len(first) >= 2:
                    return int(first[1])
    _, seq = first_fasta_record(fasta)
    return len(seq)


# -----------------------------------------------------------------------------
# Coordinate transforms and alignment-derived position maps
# -----------------------------------------------------------------------------

def rotate_pos(pos: int, p_init: int, p_len: int) -> int:
    return (int(pos) - int(p_init)) % int(p_len) + 1


def restore_human_pos(h_rot: int, offset: int, human_len: int) -> int:
    return (int(h_rot) + int(offset) - 1) % int(human_len) + 1


def parse_pairwise_alignment(aln_fa: str | Path) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    recs = read_fasta_all(aln_fa)
    if len(recs) != 2:
        die(f"Expected exactly 2 sequences in alignment, got {len(recs)}: {aln_fa}")
    (_, p_aln), (_, h_aln) = recs
    if len(p_aln) != len(h_aln):
        die("Aligned sequences have different lengths")

    cols: List[Dict[str, object]] = []
    posmap: List[Dict[str, object]] = []
    p_pos = 0
    h_pos = 0
    for i, (pc, hc) in enumerate(zip(p_aln, h_aln), start=1):
        p_is_base = pc != "-"
        h_is_base = hc != "-"
        p_curr = None
        h_curr = None
        if p_is_base:
            p_pos += 1
            p_curr = p_pos
        if h_is_base:
            h_pos += 1
            h_curr = h_pos
        if p_is_base and h_is_base:
            map_type = "base_to_base"
            if pc == hc:
                ref_relation = "match"
            elif pc in VALID and hc in VALID:
                ref_relation = "substitution"
            else:
                ref_relation = "ambiguous"
        elif p_is_base and not h_is_base:
            map_type = "human_gap"
            ref_relation = "unmapped"
        elif (not p_is_base) and h_is_base:
            map_type = "primate_gap"
            ref_relation = "NA"
        else:
            map_type = "double_gap"
            ref_relation = "NA"
        cols.append({
            "aln_col": i,
            "primate_char": pc,
            "human_char": hc,
            "primate_pos": p_curr,
            "human_pos": h_curr,
            "map_type": map_type,
            "ref_relation": ref_relation,
        })
        if p_is_base:
            posmap.append({
                "qpos": p_curr,
                "tpos": h_curr if h_is_base else "NA",
                "qref": pc,
                "href": hc if h_is_base else "NA",
                "ref_relation": ref_relation,
                "n_candidates": 1 if h_is_base else 0,
                "n_unique_tpos": 1 if h_is_base else 0,
                "continuity_flag": "NA",
                "neighbor_flag": "NA",
                "usable_for_vcf": "no",
            })
    return cols, posmap


def add_continuity_and_neighbor(posmap: List[Dict[str, object]], max_jump: int = 3, neighbor_window: int = 5) -> None:
    mapped_idx = [i for i, r in enumerate(posmap) if r["tpos"] != "NA"]
    prev = None
    for idx in mapped_idx:
        r = posmap[idx]
        if prev is None:
            r["continuity_flag"] = "start"
        else:
            dt = int(r["tpos"]) - int(prev["tpos"])
            r["continuity_flag"] = "pass" if abs(dt - 1) <= max_jump else "jump"
        prev = r
    n = len(posmap)
    for i, r in enumerate(posmap):
        if r["tpos"] == "NA":
            continue
        start = max(0, i - neighbor_window)
        end = min(n, i + neighbor_window + 1)
        nearby = [posmap[j] for j in range(start, end) if posmap[j]["tpos"] != "NA"]
        if len(nearby) < 3:
            r["neighbor_flag"] = "weak"
        else:
            jumps = sum(1 for x in nearby if x["continuity_flag"] == "jump")
            r["neighbor_flag"] = "pass" if jumps <= max(1, len(nearby) // 4) else "fail"
        if (
            r["n_unique_tpos"] == 1 and
            r["ref_relation"] == "match" and
            r["continuity_flag"] in {"pass", "start"} and
            r["neighbor_flag"] == "pass" and
            r["qref"] in VALID and
            r["href"] in VALID
        ):
            r["usable_for_vcf"] = "yes"


def write_cols(path: str | Path, cols: List[Dict[str, object]]) -> None:
    fields = ["aln_col", "primate_char", "human_char", "primate_pos", "human_pos", "map_type", "ref_relation"]
    with gzip.open(path, "wt") as out:
        out.write("\t".join(fields) + "\n")
        for r in cols:
            out.write("\t".join(str(r.get(k) if r.get(k) is not None else NA) for k in fields) + "\n")


def write_posmap(path: str | Path, posmap: List[Dict[str, object]]) -> None:
    fields = ["qpos_1based", "tpos_1based", "qref", "href", "ref_relation", "n_candidates", "n_unique_tpos", "continuity_flag", "neighbor_flag", "usable_for_vcf"]
    with gzip.open(path, "wt") as out:
        out.write("\t".join(fields) + "\n")
        for r in posmap:
            vals = [r["qpos"], r["tpos"], r["qref"], r["href"], r["ref_relation"], r["n_candidates"], r["n_unique_tpos"], r["continuity_flag"], r["neighbor_flag"], r["usable_for_vcf"]]
            out.write("\t".join(str(v) for v in vals) + "\n")


def build_map_dict(posmap_path: str | Path) -> Dict[int, Tuple[int, str, str, str]]:
    # q_rot_pos -> (human_rot_pos, usable_for_vcf, qref, href)
    mp: Dict[int, Tuple[int, str, str, str]] = {}
    with gzip.open(posmap_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[idx["tpos_1based"]] == "NA":
                continue
            mp[int(p[idx["qpos_1based"]])] = (
                int(p[idx["tpos_1based"]]),
                p[idx["usable_for_vcf"]],
                p[idx["qref"]],
                p[idx["href"]],
            )
    return mp


def write_mapping_summary(path: str | Path, sample: str, posmap: List[Dict[str, object]]) -> None:
    total = len(posmap)
    counts = Counter()
    prev_t = None
    continuous_blocks = 0
    current_block = 0
    max_block = 0
    for r in posmap:
        if r["tpos"] == "NA":
            counts["unmapped"] += 1
            if current_block > 0:
                continuous_blocks += 1
                max_block = max(max_block, current_block)
                current_block = 0
            continue
        counts["mapped"] += 1
        counts["unique"] += 1 if r["n_unique_tpos"] == 1 else 0
        counts["ambiguous"] += 0 if r["n_unique_tpos"] == 1 else 1
        counts[f"ref_{r['ref_relation']}"] += 1
        counts[f"cont_{r['continuity_flag']}"] += 1
        counts[f"neighbor_{r['neighbor_flag']}"] += 1
        if r["usable_for_vcf"] == "yes":
            counts["usable_for_vcf"] += 1
        t = int(r["tpos"])
        if prev_t is None:
            current_block = 1
        else:
            dt = t - prev_t
            if abs(dt - 1) <= 3:
                current_block += 1
            else:
                continuous_blocks += 1
                max_block = max(max_block, current_block)
                current_block = 1
        prev_t = t
    if current_block > 0:
        continuous_blocks += 1
        max_block = max(max_block, current_block)
    metrics = {
        "sample": sample,
        "total_positions": total,
        "mapped_positions": counts["mapped"],
        "unique_positions": counts["unique"],
        "ambiguous_positions": counts["ambiguous"],
        "unmapped_positions": counts["unmapped"],
        "usable_for_vcf_positions": counts["usable_for_vcf"],
        "ref_match_positions": counts["ref_match"],
        "ref_substitution_positions": counts["ref_substitution"],
        "ref_ambiguous_positions": counts["ref_ambiguous"],
        "continuity_start_positions": counts["cont_start"],
        "continuity_pass_positions": counts["cont_pass"],
        "continuity_jump_positions": counts["cont_jump"],
        "neighbor_pass_positions": counts["neighbor_pass"],
        "neighbor_fail_positions": counts["neighbor_fail"],
        "neighbor_weak_positions": counts["neighbor_weak"],
        "continuous_blocks": continuous_blocks,
        "max_continuous_block": max_block,
        "mapped_fraction": safe_div(counts["mapped"], total),
        "unique_fraction": safe_div(counts["unique"], total),
        "usable_for_vcf_fraction": safe_div(counts["usable_for_vcf"], total),
    }
    with open(path, "w") as out:
        out.write("metric\tvalue\n")
        for k, v in metrics.items():
            out.write(f"{k}\t{v}\n")


def append_summary(path: str | Path, stats: Counter, prefix: str) -> None:
    with open(path, "a") as out:
        for k, v in sorted(stats.items()):
            out.write(f"{prefix}_{k}\t{v}\n")


# -----------------------------------------------------------------------------
# VCF helpers
# -----------------------------------------------------------------------------

def parse_info(info_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info_str in {"", "."}:
        return out
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "True"
    return out


def format_info(info: Dict[str, object]) -> str:
    if not info:
        return "."
    items = []
    for k, v in info.items():
        if v is True:
            items.append(k)
        elif v is None:
            continue
        else:
            items.append(f"{k}={sanitize_info_value(v)}")
    return ";".join(items) if items else "."


def add_info_to_parts(parts: List[str], extra: Dict[str, object]) -> None:
    info = parse_info(parts[7] if len(parts) > 7 else ".")
    for k, v in extra.items():
        if v is not None:
            info[k] = v
    parts[7] = format_info(info)


def recode_gt_swap(sample_field: str) -> str:
    if sample_field in (".", "./.", ".|."):
        return sample_field
    sep = "|" if "|" in sample_field else "/"
    parts = sample_field.split(":")
    gt = parts[0]
    alleles = gt.replace("|", "/").split("/")
    new = []
    for a in alleles:
        if a == "0":
            new.append("1")
        elif a == "1":
            new.append("0")
        else:
            new.append(a)
    parts[0] = sep.join(new)
    return ":".join(parts)


def header_ids(lines: List[str], kind: str = "INFO") -> set:
    ids = set()
    prefix = f"##{kind}=<ID="
    for line in lines:
        if line.startswith(prefix):
            rest = line[len(prefix):]
            ids.add(rest.split(",", 1)[0].split(">", 1)[0])
    return ids


def write_vcf_header(fin, fout, additional_info_lines: Sequence[str], replace_chrom: Optional[str] = None) -> None:
    existing_header: List[str] = []
    added = False
    existing_info_ids = set()
    for raw in fin:
        if not raw.startswith("#"):
            # Caller must handle this record separately; not used by this helper.
            raise RuntimeError("write_vcf_header consumed a non-header line unexpectedly")
        line = raw.rstrip("\n")
        existing_header.append(line)
        if line.startswith("##INFO=<ID="):
            rest = line[len("##INFO=<ID="):]
            existing_info_ids.add(rest.split(",", 1)[0].split(">", 1)[0])
        if line.startswith("#CHROM"):
            for h in additional_info_lines:
                m = re.match(r"##INFO=<ID=([^,>]+)", h)
                if m and m.group(1) in existing_info_ids:
                    continue
                fout.write(h.rstrip("\n") + "\n")
            if replace_chrom:
                # Keep original #CHROM header unchanged. Variant CHROM values are replaced in records.
                pass
            fout.write(line + "\n")
            added = True
            break
    if not added:
        for h in additional_info_lines:
            fout.write(h.rstrip("\n") + "\n")


def stream_vcf_with_header(path: str | Path) -> Tuple[List[str], Iterator[str]]:
    # Not used in performance-critical path; helpful for small utilities.
    fh = open_text(path, "rt")
    headers = []
    def it():
        nonlocal headers
        with fh:
            for raw in fh:
                if raw.startswith("#"):
                    headers.append(raw.rstrip("\n"))
                    continue
                yield raw
    return headers, it()


def find_first_by_accession(input_dir: str | Path, accession: str, patterns: Sequence[str]) -> Optional[str]:
    hits: List[str] = []
    for pat in patterns:
        hits.extend(glob.glob(os.path.join(str(input_dir), pat.format(acc=accession))))
    hits = sorted(set(hits))
    return hits[0] if hits else None


# -----------------------------------------------------------------------------
# Alignment runner and position map stage
# -----------------------------------------------------------------------------

def run_alignment(cfg: PipeConfig, sample: str, fasta: str, aln_fasta: str, work_dir: str) -> None:
    ensure_dir(work_dir)
    pair_fasta = str(Path(work_dir) / f"{sample}.pair.fa")
    _, qseq = first_fasta_record(fasta)
    _, hseq = first_fasta_record(cfg.path_get("align_human_fasta"))
    with open(pair_fasta, "w") as out:
        out.write(f">primate\n{qseq}\n>human\n{hseq}\n")
    aligner = cfg.setting("aligner", "mafft")
    aligner_args = cfg.setting("aligner_args", "")
    if aligner == "none":
        if not Path(aln_fasta).exists():
            die(f"aligner=none but alignment does not exist: {aln_fasta}")
        return
    if shutil.which(aligner) is None:
        die(f"Aligner not found in PATH: {aligner}")
    if aligner == "mafft":
        cmd = ["mafft"] + aligner_args.split() + [pair_fasta]
        log("Running: " + " ".join(cmd))
        with open(aln_fasta, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)
    elif aligner == "muscle":
        cmd = ["muscle", "-align", pair_fasta, "-output", aln_fasta] + aligner_args.split()
        log("Running: " + " ".join(cmd))
        subprocess.run(cmd, check=True)
    elif aligner == "prank":
        prefix = str(Path(work_dir) / f"{sample}.prank")
        cmd = ["prank", f"-d={pair_fasta}", f"-o={prefix}", "-f=fasta"] + aligner_args.split()
        log("Running: " + " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copyfile(prefix + ".best.fas", aln_fasta)
    else:
        die(f"Unsupported aligner: {aligner}")


def build_posmap_stage(cfg: PipeConfig, sample: str, aln_fasta: str, dirs: Dict[str, Path]) -> Tuple[str, str, str]:
    cols, posmap = parse_pairwise_alignment(aln_fasta)
    max_jump = cfg.setting_int("max_jump", 3)
    neighbor_window = cfg.setting_int("neighbor_window", 5)
    add_continuity_and_neighbor(posmap, max_jump=max_jump, neighbor_window=neighbor_window)
    cols_path = str(dirs["maps"] / f"{sample}.aln_columns.tsv.gz")
    posmap_path = str(dirs["maps"] / f"{sample}.posmap.tsv.gz")
    summary_path = str(dirs["reports"] / f"{sample}.mapping.summary.tsv")
    write_cols(cols_path, cols)
    write_posmap(posmap_path, posmap)
    write_mapping_summary(summary_path, sample, posmap)
    return cols_path, posmap_path, summary_path


# -----------------------------------------------------------------------------
# Liftover-only VCF/COV stage
# -----------------------------------------------------------------------------

LIFTOVER_INFO_LINES = [
    '##INFO=<ID=MTLIFT_SPECIES,Number=1,Type=String,Description="Species/sample key used by mt liftover pipeline">',
    '##INFO=<ID=MTLIFT_ORIG_CHROM,Number=1,Type=String,Description="Original species chromosome/contig before liftover">',
    '##INFO=<ID=MTLIFT_ORIG_POS,Number=1,Type=Integer,Description="Original species 1-based position before liftover">',
    '##INFO=<ID=MTLIFT_ORIG_REF,Number=1,Type=String,Description="Original species REF before human REF harmonization">',
    '##INFO=<ID=MTLIFT_ORIG_ALT,Number=1,Type=String,Description="Original species ALT before human REF harmonization">',
    '##INFO=<ID=MTLIFT_ORIG_ROT_POS,Number=1,Type=Integer,Description="Original species position after species rotation">',
    '##INFO=<ID=MTLIFT_HUMAN_ROT_POS,Number=1,Type=Integer,Description="Human rotated position from pairwise alignment map">',
    '##INFO=<ID=MTLIFT_HUMAN_POS,Number=1,Type=Integer,Description="Final unrotated human position written to VCF POS">',
    '##INFO=<ID=MTLIFT_USABLE,Number=1,Type=String,Description="Whether alignment map position passed strict usable_for_vcf criteria">',
    '##INFO=<ID=MTLIFT_QREF,Number=1,Type=String,Description="Query/reference base in the rotated species alignment at mapped position">',
    '##INFO=<ID=MTLIFT_HREF_ALN,Number=1,Type=String,Description="Human base in the rotated alignment at mapped position">',
    '##INFO=<ID=MTLIFT_HUMAN_REF,Number=1,Type=String,Description="Human final reference base at lifted position">',
    '##INFO=<ID=MTLIFT_REFALT_SWAPPED,Number=1,Type=Integer,Description="1 if REF/ALT were swapped to match human reference and genotypes recoded">',
]


def make_debug_writer(debug_dir: Path, sample: str, stem: str, suffix: str, fields: Sequence[str]):
    ensure_dir(debug_dir / sample)
    path = debug_dir / sample / f"{stem}.{suffix}.tsv"
    fh = open(path, "w", newline="")
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    return fh, writer


def liftover_vcf_file(
    sample: str,
    species_key: str,
    vcf_path: str,
    out_path: str,
    drop_path: str,
    map_dict: Dict[int, Tuple[int, str, str, str]],
    sample_to_ref: Dict[str, str],
    ref_rotate_info: Dict[str, Tuple[int, int]],
    human_ref_seq: str,
    human_len: int,
    human_offset: int,
    target_chrom: str,
    debug_dir: Path,
) -> Counter:
    if sample not in sample_to_ref:
        die(f"sample={sample} not found in sample_ref_file")
    rotate_key = sample_to_ref[sample]
    if rotate_key not in ref_rotate_info:
        die(f"rotate key not found in rotate_pos_file: {rotate_key}")
    p_len, p_init = ref_rotate_info[rotate_key]
    stats = Counter()
    stem = strip_vcf_suffix(vcf_path)
    fields = [
        "sample", "species_key", "source_vcf", "reason", "original_chrom", "original_pos",
        "original_ref", "original_alt", "primate_rot_pos", "human_rot_pos", "human_final_pos",
        "usable_for_vcf", "qref", "href_aln", "human_ref_base", "written_ref", "written_alt", "filter", "info",
    ]
    fh_written, wr_written = make_debug_writer(debug_dir, sample, stem, "liftover_written", fields)
    fh_drop, wr_drop = make_debug_writer(debug_dir, sample, stem, "liftover_dropped", fields)
    try:
        with open_text(vcf_path, "rt") as fin, open_text(out_path, "wt") as fout, open_text(drop_path, "wt") as fdrop:
            # Header handling: stream manually so we do not consume first record.
            for line in fin:
                if line.startswith("##"):
                    fout.write(line)
                    fdrop.write(line)
                    continue
                if line.startswith("#CHROM"):
                    existing_ids = set()
                    # We cannot easily inspect all previous IDs here, so always append. Duplicate IDs are unlikely if raw input.
                    for h in LIFTOVER_INFO_LINES:
                        fout.write(h + "\n")
                        fdrop.write(h + "\n")
                    fout.write(line)
                    fdrop.write(line)
                    break
            for line in fin:
                if not line.strip() or line.startswith("#"):
                    continue
                stats["input_records"] += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    stats["drop_malformed"] += 1
                    fdrop.write(line)
                    continue
                orig_chrom = parts[0]
                orig_pos = int(parts[1])
                orig_ref = parts[3].upper()
                orig_alt = parts[4].upper()
                alts = orig_alt.split(",")
                p_rot = rotate_pos(orig_pos, p_init, p_len)
                hit = map_dict.get(p_rot)
                base_row = {
                    "sample": sample, "species_key": species_key, "source_vcf": os.path.basename(vcf_path),
                    "original_chrom": orig_chrom, "original_pos": orig_pos,
                    "original_ref": orig_ref, "original_alt": orig_alt,
                    "primate_rot_pos": p_rot, "filter": parts[6], "info": parts[7],
                }
                if hit is None:
                    stats["drop_unmapped"] += 1
                    fdrop.write(line)
                    row = dict(base_row, reason="unmapped", human_rot_pos=NA, human_final_pos=NA, usable_for_vcf=NA, qref=NA, href_aln=NA, human_ref_base=NA)
                    wr_drop.writerow(row)
                    continue
                h_rot, usable, qref, href_aln = hit
                h_final = restore_human_pos(h_rot, human_offset, human_len)
                base_row.update({"human_rot_pos": h_rot, "human_final_pos": h_final, "usable_for_vcf": usable, "qref": qref, "href_aln": href_aln})
                if not (len(orig_ref) == 1 and len(alts) == 1 and len(alts[0]) == 1 and orig_ref in VALID and alts[0] in VALID):
                    stats["drop_non_snv_or_multiallelic"] += 1
                    fdrop.write(line)
                    row = dict(base_row, reason="non_snv_or_multiallelic", human_ref_base=NA)
                    wr_drop.writerow(row)
                    continue
                href = human_ref_seq[h_final - 1].upper()
                base_row["human_ref_base"] = href
                if href not in VALID:
                    stats["drop_human_non_acgt"] += 1
                    fdrop.write(line)
                    wr_drop.writerow(dict(base_row, reason="human_non_acgt"))
                    continue
                parts[0] = target_chrom
                parts[1] = str(h_final)
                written_ref = orig_ref
                written_alt = alts[0]
                swapped = 0
                if href == orig_ref:
                    pass
                elif href == alts[0]:
                    parts[3], parts[4] = href, orig_ref
                    written_ref = href
                    written_alt = orig_ref
                    swapped = 1
                    for i in range(9, len(parts)):
                        parts[i] = recode_gt_swap(parts[i])
                    stats["snv_ref_alt_swapped"] += 1
                else:
                    stats["drop_ref_mismatch"] += 1
                    fdrop.write(line)
                    wr_drop.writerow(dict(base_row, reason="ref_mismatch", written_ref=NA, written_alt=NA))
                    continue
                add_info_to_parts(parts, {
                    "MTLIFT_SPECIES": species_key,
                    "MTLIFT_ORIG_CHROM": orig_chrom,
                    "MTLIFT_ORIG_POS": orig_pos,
                    "MTLIFT_ORIG_REF": orig_ref,
                    "MTLIFT_ORIG_ALT": orig_alt,
                    "MTLIFT_ORIG_ROT_POS": p_rot,
                    "MTLIFT_HUMAN_ROT_POS": h_rot,
                    "MTLIFT_HUMAN_POS": h_final,
                    "MTLIFT_USABLE": usable,
                    "MTLIFT_QREF": qref,
                    "MTLIFT_HREF_ALN": href_aln,
                    "MTLIFT_HUMAN_REF": href,
                    "MTLIFT_REFALT_SWAPPED": swapped,
                })
                fout.write("\t".join(parts) + "\n")
                stats["written_records"] += 1
                if usable == "yes":
                    stats["written_usable"] += 1
                wr_written.writerow(dict(base_row, reason="written", written_ref=written_ref, written_alt=written_alt))
    finally:
        fh_written.close()
        fh_drop.close()
    return stats


def liftover_cov_file(
    sample: str,
    cov_path: str,
    out_path: str,
    map_dict: Dict[int, Tuple[int, str, str, str]],
    sample_to_ref: Dict[str, str],
    ref_rotate_info: Dict[str, Tuple[int, int]],
    human_len: int,
    human_offset: int,
) -> Counter:
    if sample not in sample_to_ref:
        die(f"sample={sample} not found in sample_ref_file")
    rotate_key = sample_to_ref[sample]
    if rotate_key not in ref_rotate_info:
        die(f"rotate key not found in rotate_pos_file: {rotate_key}")
    p_len, p_init = ref_rotate_info[rotate_key]
    stats = Counter()
    with open_text(cov_path, "rt") as fin, open_text(out_path, "wt") as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            stats["input_rows"] += 1
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                stats["drop_malformed"] += 1
                continue
            try:
                pos = int(p[1])
            except ValueError:
                stats["drop_bad_pos"] += 1
                continue
            p_rot = rotate_pos(pos, p_init, p_len)
            hit = map_dict.get(p_rot)
            if hit is None:
                stats["drop_unmapped"] += 1
                continue
            h_rot, usable, _qref, _href = hit
            h_final = restore_human_pos(h_rot, human_offset, human_len)
            p[1] = str(h_final)
            fout.write("\t".join(p) + "\n")
            stats["written_rows"] += 1
            if usable == "yes":
                stats["written_usable"] += 1
    return stats


def liftover_stage(cfg: PipeConfig, sample: str, posmap_path: str, summary_path: str, dirs: Dict[str, Path]) -> List[str]:
    map_dict = build_map_dict(posmap_path)
    human_ref_path = cfg.path_get("ref_human_fasta")
    recs = read_fasta_all(human_ref_path)
    if len(recs) != 1:
        warn(f"Human REF FASTA has {len(recs)} records; using the first")
    human_ref_seq = recs[0][1]
    sample_to_ref = load_sample_ref_map(cfg.path_get("sample_ref_file"))
    species_to_accessions = load_species_to_accessions(cfg.path_get("sample_ref_file"))
    ref_rotate_info = load_rotate_info(cfg.path_get("rotate_pos_file"))
    human_len = cfg.setting_int("human_len", 16569)
    human_offset = cfg.setting_int("human_restore_offset", 1325)
    target_chrom = cfg.setting("target_chrom", "chrM")
    raw_vcfs: List[str] = []
    accession_list = species_to_accessions.get(sample, [])
    if cfg.setting_bool("run_vcf_liftover", True):
        sample_out_dir = ensure_dir(dirs["vcf_raw"] / sample)
        vcf_input_dir = cfg.path_get("vcf_input_dir")
        found_n = 0
        agg = Counter(samples_expected=len(accession_list))
        for accession in accession_list:
            vcf_in = find_first_by_accession(vcf_input_dir, accession, ["{acc}*.vcf", "{acc}*.vcf.gz"])
            if vcf_in is None:
                warn(f"No VCF found for sample={sample}, accession={accession}")
                continue
            found_n += 1
            stem = strip_vcf_suffix(vcf_in)
            out_vcf = str(sample_out_dir / f"{stem}.lifted.raw.vcf")
            drop_vcf = str(sample_out_dir / f"{stem}.lifted.raw.dropped.vcf")
            log(f"Liftover VCF: {vcf_in} -> {out_vcf}")
            stats = liftover_vcf_file(sample, sample, vcf_in, out_vcf, drop_vcf, map_dict, sample_to_ref, ref_rotate_info, human_ref_seq, human_len, human_offset, target_chrom, dirs["debug"])
            agg.update(stats)
            raw_vcfs.append(out_vcf)
        agg["samples_found"] = found_n
        agg["file_found"] = 1 if found_n > 0 else 0
        append_summary(summary_path, agg, "vcf_liftover")
    if cfg.setting_bool("run_cov_liftover", False):
        sample_out_dir = ensure_dir(dirs["cov"] / sample)
        cov_input_dir = cfg.path_get("cov_input_dir")
        found_n = 0
        agg = Counter(samples_expected=len(accession_list))
        for accession in accession_list:
            cov_in = find_first_by_accession(cov_input_dir, accession, ["{acc}*.tsv", "{acc}*.tsv.gz"])
            if cov_in is None:
                warn(f"No COV found for sample={sample}, accession={accession}")
                continue
            found_n += 1
            stem = Path(cov_in).name
            if stem.endswith(".gz"):
                stem = stem[:-3]
            out_cov = str(sample_out_dir / f"{stem}.lifted.tsv")
            stats = liftover_cov_file(sample, cov_in, out_cov, map_dict, sample_to_ref, ref_rotate_info, human_len, human_offset)
            agg.update(stats)
        agg["samples_found"] = found_n
        agg["file_found"] = 1 if found_n > 0 else 0
        append_summary(summary_path, agg, "cov_liftover")
    return raw_vcfs


# -----------------------------------------------------------------------------
# Codon match annotation-only stage
# -----------------------------------------------------------------------------

def harmonize_gene_name(g: str) -> str:
    g = g.strip()
    mapping = {
        "CYB": "MT-CYB", "CYTB": "MT-CYB",
        "ATP6": "MT-ATP6", "ATP8": "MT-ATP8",
        "CO1": "MT-CO1", "COI": "MT-CO1", "COX1": "MT-CO1",
        "CO2": "MT-CO2", "COII": "MT-CO2", "COX2": "MT-CO2",
        "CO3": "MT-CO3", "COIII": "MT-CO3", "COX3": "MT-CO3",
        "ND1": "MT-ND1", "ND2": "MT-ND2", "ND3": "MT-ND3",
        "ND4": "MT-ND4", "ND4L": "MT-ND4L", "ND5": "MT-ND5", "ND6": "MT-ND6",
    }
    return mapping.get(g, g)


def load_all_primate_position_codon_table(path: str | Path) -> Dict[Tuple[str, int], List[Dict[str, object]]]:
    out: Dict[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}
        need = ["pos", "gene", "strand", "codon_index", "codon_pos_in_triplet", "codon_seq", "codon_pos1_genomic", "codon_pos2_genomic", "codon_pos3_genomic"]
        for req in need:
            if req not in idx:
                die(f"{path} missing required column: {req}")
        id_cols = [c for c in ["species_key", "species", "file_name", "seq_name", "accession"] if c in idx]
        if not id_cols:
            die(f"{path} needs one identifier column: species_key/species/file_name/seq_name/accession")
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[idx["pos"]])
            row = {
                "gene": harmonize_gene_name(p[idx["gene"]]),
                "strand": p[idx["strand"]],
                "codon_index": int(p[idx["codon_index"]]),
                "codon_pos_in_triplet": int(p[idx["codon_pos_in_triplet"]]),
                "codon_seq": p[idx["codon_seq"]].upper(),
                "codon_pos1_genomic": int(p[idx["codon_pos1_genomic"]]),
                "codon_pos2_genomic": int(p[idx["codon_pos2_genomic"]]),
                "codon_pos3_genomic": int(p[idx["codon_pos3_genomic"]]),
            }
            aliases = set()
            for col in id_cols:
                val = p[idx[col]]
                if val and val != "NA":
                    aliases.add(val)
                    if col in {"species", "file_name"}:
                        aliases.add(val.replace(" ", "_"))
                        aliases.add(os.path.splitext(val)[0])
            for a in aliases:
                out[(a, pos)].append(dict(row))
    return out


def load_human_codon_lookup(path: str | Path) -> Dict[int, List[Dict[str, object]]]:
    out: Dict[int, List[Dict[str, object]]] = defaultdict(list)
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}
        need = ["pos", "gene", "strand", "codon_index", "codon_pos_in_triplet", "codon_seq", "codon_pos1_genomic", "codon_pos2_genomic", "codon_pos3_genomic"]
        for req in need:
            if req not in idx:
                die(f"{path} missing required column: {req}")
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            pos = int(p[idx["pos"]])
            out[pos].append({
                "gene": harmonize_gene_name(p[idx["gene"]]),
                "strand": p[idx["strand"]],
                "codon_index": int(p[idx["codon_index"]]),
                "codon_pos_in_triplet": int(p[idx["codon_pos_in_triplet"]]),
                "codon_seq": p[idx["codon_seq"]].upper(),
                "codon_pos1_genomic": int(p[idx["codon_pos1_genomic"]]),
                "codon_pos2_genomic": int(p[idx["codon_pos2_genomic"]]),
                "codon_pos3_genomic": int(p[idx["codon_pos3_genomic"]]),
            })
    return out


def mutate_codon(codon: str, pos_in_triplet: int, alt_base: str) -> str:
    codon_list = list(codon)
    codon_list[int(pos_in_triplet) - 1] = alt_base.upper()
    return "".join(codon_list)


def choose_human_codon_candidates(primate_row: Dict[str, object], human_rows: List[Dict[str, object]], strict_phase_match: bool = False) -> List[Dict[str, object]]:
    """Return human codon rows that are eligible for codon-sequence comparison.

    In strict phase mode, a candidate must be in the same mitochondrial gene AND
    the lifted human position must occupy the same codon phase (1/2/3) as the
    source/species variant. This prevents a codon-sequence match from passing
    when the variant is, for example, base 1 of a primate codon but base 2 of
    the corresponding human codon.

    In non-strict mode, we keep the older permissive behavior: prefer same-gene
    candidates, otherwise fall back to all human coding rows at the lifted site.
    """
    if not human_rows:
        return []
    gene = harmonize_gene_name(str(primate_row["gene"]))
    same_gene = [r for r in human_rows if harmonize_gene_name(str(r["gene"])) == gene]
    if strict_phase_match:
        phase = int(primate_row["codon_pos_in_triplet"])
        return [r for r in same_gene if int(r["codon_pos_in_triplet"]) == phase]
    if same_gene:
        return same_gene
    return human_rows


def summarize_human_codon_compatibility(
    primate_rows: List[Dict[str, object]],
    human_rows: List[Dict[str, object]],
) -> Tuple[bool, bool]:
    """Return whether any human row matches source gene and source codon phase."""
    any_same_gene = False
    any_same_phase_within_same_gene = False
    for prow in primate_rows:
        pgene = harmonize_gene_name(str(prow["gene"]))
        pphase = int(prow["codon_pos_in_triplet"])
        for hrow in human_rows:
            if harmonize_gene_name(str(hrow["gene"])) != pgene:
                continue
            any_same_gene = True
            if int(hrow["codon_pos_in_triplet"]) == pphase:
                any_same_phase_within_same_gene = True
    return any_same_gene, any_same_phase_within_same_gene


CODON_INFO_LINES = [
    '##INFO=<ID=MTCODON_STATUS,Number=1,Type=String,Description="Codon-match annotation status: PASS, SKIPPED_NONCODING, NO_HUMAN_CODON, GENE_MISMATCH, PHASE_MISMATCH, MISMATCH, or MISSING_COORD">',
    '##INFO=<ID=MTCODON_MATCH,Number=1,Type=String,Description="How codon matched human codon after required gene/phase checks: ref_codon, alt_codon, or NA">',
    '##INFO=<ID=MTCODON_STRICT_PHASE,Number=1,Type=Integer,Description="Whether strict codon phase matching was required: 1=yes, 0=no">',
    '##INFO=<ID=MTCODON_GENE_MATCH,Number=1,Type=String,Description="Whether primate/source and human codon annotations are in the same mitochondrial gene: yes, no, or NA">',
    '##INFO=<ID=MTCODON_PHASE_MATCH,Number=1,Type=String,Description="Whether primate/source variant codon position in triplet equals the human lifted position codon position in triplet: yes, no, or NA">',
    '##INFO=<ID=MTCODON_PRIMATE_GENE,Number=1,Type=String,Description="Primate/source mitochondrial coding gene">',
    '##INFO=<ID=MTCODON_PRIMATE_CODON,Number=1,Type=String,Description="Primate/source codon sequence">',
    '##INFO=<ID=MTCODON_PRIMATE_PHASE,Number=1,Type=Integer,Description="Primate/source position in codon triplet, 1-3">',
    '##INFO=<ID=MTCODON_HUMAN_GENE,Number=1,Type=String,Description="Human mitochondrial coding gene at lifted position">',
    '##INFO=<ID=MTCODON_HUMAN_CODON,Number=1,Type=String,Description="Human codon sequence at lifted position">',
    '##INFO=<ID=MTCODON_HUMAN_PHASE,Number=1,Type=Integer,Description="Human position in codon triplet, 1-3">',
]


def representative_human_codon_row(
    primate_rows: List[Dict[str, object]],
    human_rows: List[Dict[str, object]],
) -> Dict[str, object]:
    """Pick the most informative human row for diagnostic INFO on failures."""
    for prow in primate_rows:
        for hrow in human_rows:
            if harmonize_gene_name(str(prow["gene"])) == harmonize_gene_name(str(hrow["gene"])) and int(prow["codon_pos_in_triplet"]) == int(hrow["codon_pos_in_triplet"]):
                return hrow
    for prow in primate_rows:
        for hrow in human_rows:
            if harmonize_gene_name(str(prow["gene"])) == harmonize_gene_name(str(hrow["gene"])):
                return hrow
    return human_rows[0]


def yesno(x: bool) -> str:
    return "yes" if x else "no"


def codon_match_for_record(
    parts: List[str],
    species_key: str,
    primate_lookup: Dict[Tuple[str, int], List[Dict[str, object]]],
    human_lookup: Dict[int, List[Dict[str, object]]],
    strict_phase_match: bool,
) -> Dict[str, object]:
    info = parse_info(parts[7])
    orig_pos = parse_optional_int(info.get("MTLIFT_ORIG_POS"))
    human_pos = parse_optional_int(info.get("MTLIFT_HUMAN_POS")) or int(parts[1])
    orig_alt = info.get("MTLIFT_ORIG_ALT", parts[4]).split(",")[0].upper()
    base_ann = {
        "MTCODON_STRICT_PHASE": int(bool(strict_phase_match)),
        "MTCODON_GENE_MATCH": NA,
        "MTCODON_PHASE_MATCH": NA,
    }
    if orig_pos is None or human_pos is None:
        return {**base_ann, "MTCODON_STATUS": "MISSING_COORD", "MTCODON_MATCH": NA}

    primate_rows = primate_lookup.get((species_key, orig_pos), [])
    if not primate_rows:
        return {**base_ann, "MTCODON_STATUS": "SKIPPED_NONCODING", "MTCODON_MATCH": NA}

    human_rows = human_lookup.get(human_pos, [])
    if not human_rows:
        return {
            **base_ann,
            "MTCODON_STATUS": "NO_HUMAN_CODON",
            "MTCODON_MATCH": NA,
            "MTCODON_PRIMATE_GENE": primate_rows[0]["gene"],
            "MTCODON_PRIMATE_CODON": primate_rows[0]["codon_seq"],
            "MTCODON_PRIMATE_PHASE": primate_rows[0]["codon_pos_in_triplet"],
        }

    any_same_gene, any_same_phase_within_same_gene = summarize_human_codon_compatibility(primate_rows, human_rows)

    for prow in primate_rows:
        candidates = choose_human_codon_candidates(prow, human_rows, strict_phase_match=strict_phase_match)
        if not candidates:
            continue
        primate_ref_codon = str(prow["codon_seq"]).upper()
        primate_alt_codon = mutate_codon(primate_ref_codon, int(prow["codon_pos_in_triplet"]), orig_alt)
        for hrow in candidates:
            human_codon = str(hrow["codon_seq"]).upper()
            matched_by = None
            if primate_ref_codon == human_codon:
                matched_by = "ref_codon"
            elif primate_alt_codon == human_codon:
                matched_by = "alt_codon"
            if matched_by:
                gene_match = harmonize_gene_name(str(prow["gene"])) == harmonize_gene_name(str(hrow["gene"]))
                phase_match = int(prow["codon_pos_in_triplet"]) == int(hrow["codon_pos_in_triplet"])
                # In strict mode, PASS is allowed only for same-gene + same-phase candidates.
                # This guard is intentionally redundant with choose_human_codon_candidates() so
                # future edits cannot accidentally reintroduce permissive fallback behavior.
                if strict_phase_match and not (gene_match and phase_match):
                    continue
                return {
                    **base_ann,
                    "MTCODON_STATUS": "PASS",
                    "MTCODON_MATCH": matched_by,
                    "MTCODON_GENE_MATCH": yesno(gene_match),
                    "MTCODON_PHASE_MATCH": yesno(phase_match),
                    "MTCODON_PRIMATE_GENE": prow["gene"],
                    "MTCODON_PRIMATE_CODON": prow["codon_seq"],
                    "MTCODON_PRIMATE_PHASE": prow["codon_pos_in_triplet"],
                    "MTCODON_HUMAN_GENE": hrow["gene"],
                    "MTCODON_HUMAN_CODON": hrow["codon_seq"],
                    "MTCODON_HUMAN_PHASE": hrow["codon_pos_in_triplet"],
                }

    diagnostic_hrow = representative_human_codon_row(primate_rows, human_rows)
    status = "MISMATCH"
    if strict_phase_match and not any_same_gene:
        status = "GENE_MISMATCH"
    elif strict_phase_match and not any_same_phase_within_same_gene:
        status = "PHASE_MISMATCH"
    gene_match = any_same_gene
    phase_match = any_same_phase_within_same_gene
    return {
        **base_ann,
        "MTCODON_STATUS": status,
        "MTCODON_MATCH": NA,
        "MTCODON_GENE_MATCH": yesno(gene_match),
        "MTCODON_PHASE_MATCH": yesno(phase_match),
        "MTCODON_PRIMATE_GENE": primate_rows[0]["gene"],
        "MTCODON_PRIMATE_CODON": primate_rows[0]["codon_seq"],
        "MTCODON_PRIMATE_PHASE": primate_rows[0]["codon_pos_in_triplet"],
        "MTCODON_HUMAN_GENE": diagnostic_hrow["gene"],
        "MTCODON_HUMAN_CODON": diagnostic_hrow["codon_seq"],
        "MTCODON_HUMAN_PHASE": diagnostic_hrow["codon_pos_in_triplet"],
    }


def annotate_codon_vcf_file(
    in_vcf: str,
    out_vcf: str,
    species_key: str,
    primate_lookup: Dict[Tuple[str, int], List[Dict[str, object]]],
    human_lookup: Dict[int, List[Dict[str, object]]],
    strict_phase_match: bool,
    summary_path: Optional[str] = None,
) -> Counter:
    stats = Counter()
    with open_text(in_vcf, "rt") as fin, open_text(out_vcf, "wt") as fout:
        header_done = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                for h in CODON_INFO_LINES:
                    fout.write(h + "\n")
                fout.write(line)
                header_done = True
                break
        if not header_done:
            die(f"Malformed VCF header: {in_vcf}")
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                stats["malformed"] += 1
                continue
            ann = codon_match_for_record(parts, species_key, primate_lookup, human_lookup, strict_phase_match)
            add_info_to_parts(parts, ann)
            status = str(ann.get("MTCODON_STATUS", "UNKNOWN"))
            stats[f"status_{status}"] += 1
            stats["written_records"] += 1
            fout.write("\t".join(parts) + "\n")
    if summary_path:
        append_summary(summary_path, stats, "codon")
    return stats


def codon_stage(cfg: PipeConfig, sample: str, input_vcfs: List[str], summary_path: str, dirs: Dict[str, Path]) -> List[str]:
    primate_table = cfg.path_get("all_primate_position_codon_table")
    human_table = cfg.path_get("human_codon_table")
    primate_lookup = load_all_primate_position_codon_table(primate_table)
    human_lookup = load_human_codon_lookup(human_table)
    strict = cfg.setting_bool("strict_phase_match", True)
    out_dir = ensure_dir(dirs["vcf_codon"] / sample)
    outputs: List[str] = []
    for in_vcf in input_vcfs:
        stem = strip_vcf_suffix(in_vcf)
        # Remove .lifted.raw suffix if present for cleaner names.
        if stem.endswith(".lifted.raw"):
            stem = stem[: -len(".lifted.raw")]
        out_vcf = str(out_dir / f"{stem}.lifted.codon.vcf")
        log(f"Codon annotate: {in_vcf} -> {out_vcf}")
        annotate_codon_vcf_file(in_vcf, out_vcf, sample, primate_lookup, human_lookup, strict, summary_path)
        outputs.append(out_vcf)
    return outputs


# -----------------------------------------------------------------------------
# tRNA scan, index, and match annotation-only stage
# -----------------------------------------------------------------------------

OUT_HEADER_RE = re.compile(r"^\s*(Sequence|Name)\s+", re.I)
DASH_RE = re.compile(r"^-{3,}")
SS_HEADER_RE = re.compile(r"^(?P<trna_id>\S+)\s+\((?P<begin>\d+)\s*-\s*(?P<end>\d+)\)\s+Length:\s+(?P<length>\d+)\s+bp", re.I)
SS_TYPE_RE = re.compile(r"Type:\s+(?P<aa>\S+)\s+Anticodon:\s+(?P<anticodon>\S+)\s+at\s+(?P<anti_begin>\d+)\s*-\s*(?P<anti_end>\d+)", re.I)


@dataclass
class TRNARecord:
    trna_id: str
    chrom: str
    trna_num: str
    begin_raw: int
    end_raw: int
    start: int
    end: int
    strand: str
    aa: str = NA
    anticodon: str = NA
    intron_begin: Optional[int] = None
    intron_end: Optional[int] = None
    score: Optional[float] = None
    seq: str = ""
    struct: str = ""
    pairs: Dict[int, int] = field(default_factory=dict)
    elements: Dict[int, str] = field(default_factory=dict)

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end

    def local_index(self, genome_pos: int) -> Optional[int]:
        if not self.contains(genome_pos):
            return None
        if self.strand == "+":
            return genome_pos - self.begin_raw + 1
        return self.begin_raw - genome_pos + 1

    def genome_pos(self, local_idx: Optional[int]) -> Optional[int]:
        if local_idx is None or local_idx < 1:
            return None
        n = max(len(self.seq), abs(self.end_raw - self.begin_raw) + 1)
        if local_idx > n:
            return None
        if self.strand == "+":
            return self.begin_raw + local_idx - 1
        return self.begin_raw - local_idx + 1

    def base_at(self, local_idx: Optional[int]) -> str:
        if local_idx is None or local_idx < 1 or local_idx > len(self.seq):
            return NA
        return self.seq[local_idx - 1].upper().replace("T", "U")

    def struct_char_at(self, local_idx: Optional[int]) -> str:
        if local_idx is None or local_idx < 1 or local_idx > len(self.struct):
            return NA
        return self.struct[local_idx - 1]

    def struct_class_at(self, local_idx: Optional[int]) -> str:
        c = self.struct_char_at(local_idx)
        if c in {">", "<"}:
            return "stem"
        if c == ".":
            return "loop"
        return NA

    def struct_element_at(self, local_idx: Optional[int]) -> str:
        if local_idx is None:
            return NA
        return self.elements.get(local_idx, self.struct_class_at(local_idx))

    def pair_index(self, local_idx: Optional[int]) -> Optional[int]:
        if local_idx is None:
            return None
        return self.pairs.get(local_idx)


def normalize_chrom(chrom: Optional[str], mode: str) -> Optional[str]:
    if chrom is None:
        return None
    chrom = str(chrom)
    if mode == "strip_chr" and chrom.startswith("chr"):
        return chrom[3:]
    if mode == "add_chr" and not chrom.startswith("chr"):
        return "chr" + chrom
    return chrom


def parse_trnascan_out(out_path: str | Path, chrom_norm: str = "none") -> Dict[str, TRNARecord]:
    records: Dict[str, TRNARecord] = {}
    with open_text(out_path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or OUT_HEADER_RE.match(line) or DASH_RE.match(line):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            chrom = normalize_chrom(parts[0], chrom_norm) or parts[0]
            trna_num = parts[1]
            begin = parse_optional_int(parts[2])
            end = parse_optional_int(parts[3])
            if begin is None or end is None:
                continue
            aa = parts[4] if len(parts) > 4 else NA
            anticodon = parts[5] if len(parts) > 5 else NA
            intron_begin = parse_optional_int(parts[6]) if len(parts) > 6 else None
            intron_end = parse_optional_int(parts[7]) if len(parts) > 7 else None
            score = parse_optional_float(parts[8]) if len(parts) > 8 else None
            strand = "+" if begin <= end else "-"
            trna_id = f"{chrom}.trna{trna_num}"
            records[trna_id] = TRNARecord(
                trna_id=trna_id, chrom=chrom, trna_num=trna_num,
                begin_raw=begin, end_raw=end, start=min(begin, end), end=max(begin, end),
                strand=strand, aa=aa, anticodon=anticodon,
                intron_begin=intron_begin, intron_end=intron_end, score=score,
            )
    return records


def build_pairs(struct: str) -> Dict[int, int]:
    stack: List[int] = []
    pairs: Dict[int, int] = {}
    for idx, char in enumerate(struct, start=1):
        if char == ">":
            stack.append(idx)
        elif char == "<":
            if not stack:
                continue
            left = stack.pop()
            pairs[left] = idx
            pairs[idx] = left
    return pairs


def _stem_groups_from_struct(struct: str, pairs: Dict[int, int]) -> List[Dict[str, int]]:
    groups: List[Dict[str, int]] = []
    n = len(struct)
    i = 1
    while i <= n:
        if struct[i - 1] != ">":
            i += 1
            continue
        left_start = i
        while i <= n and struct[i - 1] == ">":
            i += 1
        left_end = i - 1
        right = [pairs.get(p) for p in range(left_start, left_end + 1) if pairs.get(p) is not None]
        if right:
            groups.append({"left_start": left_start, "left_end": left_end, "right_start": min(right), "right_end": max(right)})
    return groups


def infer_structural_elements(struct: str, pairs: Dict[int, int]) -> Dict[int, str]:
    n = len(struct)
    elements: Dict[int, str] = {i: ("stem" if struct[i - 1] in "><" else "loop") for i in range(1, n + 1)}
    groups = _stem_groups_from_struct(struct, pairs)
    if not groups:
        return elements
    acceptor_idx = None
    rightmost = -1
    for idx, g in enumerate(groups):
        if g["right_end"] > rightmost:
            rightmost = g["right_end"]
            acceptor_idx = idx
    labels_by_group: Dict[int, str] = {}
    if acceptor_idx is not None:
        labels_by_group[acceptor_idx] = "acceptor_stem"
    internal = [i for i in range(len(groups)) if i != acceptor_idx]
    internal.sort(key=lambda k: groups[k]["left_start"])
    if internal:
        labels_by_group[internal[0]] = "D_stem"
    if len(internal) >= 2:
        labels_by_group[internal[1]] = "anticodon_stem"
    if len(internal) == 3:
        labels_by_group[internal[2]] = "T_stem"
    elif len(internal) >= 4:
        for k in internal[2:-1]:
            labels_by_group[k] = "variable_stem"
        labels_by_group[internal[-1]] = "T_stem"
    for idx, g in enumerate(groups):
        label = labels_by_group.get(idx, "stem")
        for p in range(g["left_start"], g["left_end"] + 1):
            elements[p] = label
        for p in range(g["right_start"], g["right_end"] + 1):
            elements[p] = label
        if label != "acceptor_stem":
            loop_label = label.replace("_stem", "_loop") if label.endswith("_stem") else "loop"
            for p in range(g["left_end"] + 1, g["right_start"]):
                if 1 <= p <= n and struct[p - 1] == ".":
                    elements[p] = loop_label
    return elements


def _finalize_ss_record(current: Optional[Dict[str, object]], trnas: Dict[str, TRNARecord], chrom_norm: str) -> None:
    if not current:
        return
    raw_id = str(current.get("trna_id"))
    chrom0 = raw_id.split(".trna", 1)[0]
    chrom = normalize_chrom(chrom0, chrom_norm) or chrom0
    m_num = re.search(r"\.trna(\S+)$", raw_id)
    normalized_id = f"{chrom}.trna{m_num.group(1)}" if m_num else raw_id
    rec = trnas.get(normalized_id)
    if rec is None:
        begin = int(current.get("begin", 0))
        end = int(current.get("end", 0))
        strand = "+" if begin <= end else "-"
        rec = TRNARecord(normalized_id, chrom, m_num.group(1) if m_num else NA, begin, end, min(begin, end), max(begin, end), strand)
        trnas[normalized_id] = rec
    seq = str(current.get("seq", "")).replace(" ", "").upper().replace("T", "U")
    struct = str(current.get("struct", "")).replace(" ", "")
    rec.seq = seq
    rec.struct = struct
    rec.pairs = build_pairs(struct)
    rec.elements = infer_structural_elements(struct, rec.pairs)
    if current.get("aa"):
        rec.aa = str(current.get("aa"))
    if current.get("anticodon"):
        rec.anticodon = str(current.get("anticodon"))


def parse_trnascan_ss(ss_path: str | Path, trnas: Dict[str, TRNARecord], chrom_norm: str = "none") -> Dict[str, TRNARecord]:
    current: Optional[Dict[str, object]] = None
    with open_text(ss_path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            m = SS_HEADER_RE.match(line)
            if m:
                _finalize_ss_record(current, trnas, chrom_norm)
                current = {"trna_id": m.group("trna_id"), "begin": int(m.group("begin")), "end": int(m.group("end")), "length": int(m.group("length"))}
                rest = line[m.end():].strip()
                if rest:
                    mt = SS_TYPE_RE.search(rest)
                    if mt:
                        current["aa"] = mt.group("aa")
                        current["anticodon"] = mt.group("anticodon")
                    ms = re.search(r"Seq:\s*([A-Za-z]+)", rest)
                    if ms:
                        current["seq"] = ms.group(1)
                    mst = re.search(r"Str:\s*([<>.]+)", rest)
                    if mst:
                        current["struct"] = mst.group(1)
                continue
            if current is None:
                continue
            mt = SS_TYPE_RE.search(line)
            if mt:
                current["aa"] = mt.group("aa")
                current["anticodon"] = mt.group("anticodon")
            elif line.startswith("Seq:"):
                current["seq"] = line.split("Seq:", 1)[1].strip().replace(" ", "")
            elif line.startswith("Str:"):
                current["struct"] = line.split("Str:", 1)[1].strip().replace(" ", "")
            else:
                ms = re.search(r"Seq:\s*([A-Za-z]+)", line)
                if ms:
                    current["seq"] = ms.group(1)
                mst = re.search(r"Str:\s*([<>.]+)", line)
                if mst:
                    current["struct"] = mst.group(1)
    _finalize_ss_record(current, trnas, chrom_norm)
    return trnas


def merge_trnas(out_path: str | Path, ss_path: str | Path, chrom_norm: str = "none") -> Dict[str, TRNARecord]:
    trnas = parse_trnascan_out(out_path, chrom_norm=chrom_norm)
    trnas = parse_trnascan_ss(ss_path, trnas, chrom_norm=chrom_norm)
    missing = [r.trna_id for r in trnas.values() if not r.seq or not r.struct]
    if missing:
        warn(f"{len(missing)} tRNA records lack Seq/Str structure from {ss_path}")
    return trnas


def normalize_base(base: str) -> str:
    if base is None or base == NA:
        return NA
    b = str(base).upper().replace("T", "U")
    if len(b) != 1 or b not in "ACGUN":
        return NA
    return b


def pair_type(base1: str, base2: str) -> str:
    b1 = normalize_base(base1)
    b2 = normalize_base(base2)
    if b1 == NA or b2 == NA or "N" in {b1, b2}:
        return NA
    pair = b1 + b2
    if pair in {"AU", "UA", "GC", "CG"}:
        return "WC"
    if pair in {"GU", "UG"}:
        return "GU_wobble"
    return "non_WC"


def binary_pair_state(ptype: str) -> str:
    if ptype == NA:
        return NA
    return "WC" if ptype == "WC" else "non_WC"


def trna_index_fields() -> List[str]:
    return [
        "species_key", "chrom", "pos", "trna_id", "trna_begin", "trna_end", "strand", "aa", "anticodon", "score",
        "local_pos", "base", "struct_char", "struct_class", "struct_element",
        "paired_local_pos", "paired_genomic_pos", "paired_base", "pair_bases", "pair_type", "pair_state",
    ]


def build_trna_position_index(species_key: str, trnascan_out: str, trnascan_ss: str, output_tsv: str, chrom_norm: str = "none") -> None:
    records = merge_trnas(trnascan_out, trnascan_ss, chrom_norm=chrom_norm)
    fields = trna_index_fields()
    rows: List[Dict[str, object]] = []
    for rec in records.values():
        n = max(len(rec.seq), abs(rec.end_raw - rec.begin_raw) + 1)
        for local in range(1, n + 1):
            pos = rec.genome_pos(local)
            if pos is None:
                continue
            pair_local = rec.pair_index(local)
            pair_genome = rec.genome_pos(pair_local) if pair_local is not None else None
            base = rec.base_at(local)
            pbase = rec.base_at(pair_local)
            ptype = pair_type(base, pbase) if pair_local is not None else NA
            rows.append({
                "species_key": species_key,
                "chrom": rec.chrom,
                "pos": pos,
                "trna_id": rec.trna_id,
                "trna_begin": rec.begin_raw,
                "trna_end": rec.end_raw,
                "strand": rec.strand,
                "aa": rec.aa,
                "anticodon": rec.anticodon,
                "score": rec.score if rec.score is not None else NA,
                "local_pos": local,
                "base": base,
                "struct_char": rec.struct_char_at(local),
                "struct_class": rec.struct_class_at(local),
                "struct_element": rec.struct_element_at(local),
                "paired_local_pos": pair_local if pair_local is not None else NA,
                "paired_genomic_pos": pair_genome if pair_genome is not None else NA,
                "paired_base": pbase,
                "pair_bases": f"{normalize_base(base)}-{normalize_base(pbase)}" if pair_local is not None else NA,
                "pair_type": ptype,
                "pair_state": binary_pair_state(ptype),
            })
    with open_text(output_tsv, "wt") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in sorted(rows, key=lambda r: (str(r["chrom"]), int(r["pos"]), str(r["trna_id"]))):
            w.writerow(row)
    log(f"Wrote tRNA position index: {output_tsv}")


def load_trna_position_index(
    path: str | Path,
    chrom_norm: str = "none",
    ignore_chrom: bool = False,
) -> Dict[Tuple[str, int], List[Dict[str, str]]]:
    idx: Dict[Tuple[str, int], List[Dict[str, str]]] = defaultdict(list)
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = normalize_chrom(row.get("chrom"), chrom_norm)
            pos = parse_optional_int(row.get("pos"))
            if pos is None:
                continue
            if ignore_chrom:
                key = ("*", pos)
            else:
                if chrom is None:
                    continue
                key = (chrom, pos)
            idx[key].append(row)
    for k in idx:
        idx[k].sort(key=lambda r: r.get("trna_id", ""))
    return idx


TRNA_INFO_LINES = [
    '##INFO=<ID=MTTRNA_STATUS,Number=1,Type=String,Description="tRNA annotation status">',
    '##INFO=<ID=MTTRNA_S_ID,Number=1,Type=String,Description="Species/source tRNA ID overlapping original coordinate">',
    '##INFO=<ID=MTTRNA_H_ID,Number=1,Type=String,Description="Human tRNA ID overlapping lifted coordinate">',
    '##INFO=<ID=MTTRNA_S_LOCAL,Number=1,Type=String,Description="Species/source local tRNA position">',
    '##INFO=<ID=MTTRNA_H_LOCAL,Number=1,Type=String,Description="Human local tRNA position">',
    '##INFO=<ID=MTTRNA_S_CLASS,Number=1,Type=String,Description="Species/source structural class: stem or loop">',
    '##INFO=<ID=MTTRNA_H_CLASS,Number=1,Type=String,Description="Human structural class: stem or loop">',
    '##INFO=<ID=MTTRNA_REGION_MATCH,Number=1,Type=String,Description="yes/no/NA, whether source and human stem/loop class match">',
    '##INFO=<ID=MTTRNA_S_ELEMENT,Number=1,Type=String,Description="Species/source structural element label">',
    '##INFO=<ID=MTTRNA_H_ELEMENT,Number=1,Type=String,Description="Human structural element label">',
    '##INFO=<ID=MTTRNA_ELEMENT_MATCH,Number=1,Type=String,Description="yes/no/NA, whether source and human structural element labels match">',
    '##INFO=<ID=MTTRNA_S_PAIR_TYPE,Number=1,Type=String,Description="Species/source pair type: WC, GU_wobble, non_WC, or NA">',
    '##INFO=<ID=MTTRNA_H_PAIR_TYPE,Number=1,Type=String,Description="Human pair type: WC, GU_wobble, non_WC, or NA">',
    '##INFO=<ID=MTTRNA_PAIR_TYPE_MATCH,Number=1,Type=String,Description="yes/no/NA, whether source and human pair_type match">',
    '##INFO=<ID=MTTRNA_S_PAIR_STATE,Number=1,Type=String,Description="Species/source binary pair state: WC, non_WC, or NA">',
    '##INFO=<ID=MTTRNA_H_PAIR_STATE,Number=1,Type=String,Description="Human binary pair state: WC, non_WC, or NA">',
    '##INFO=<ID=MTTRNA_PAIR_STATE_MATCH,Number=1,Type=String,Description="yes/no/NA, whether binary pair state matches">',
    '##INFO=<ID=MTTRNA_S_PAIR_LOCAL,Number=1,Type=String,Description="Species/source paired local tRNA position">',
    '##INFO=<ID=MTTRNA_H_PAIR_LOCAL,Number=1,Type=String,Description="Human paired local tRNA position">',
    '##INFO=<ID=MTTRNA_PAIR_LOCAL_MATCH,Number=1,Type=String,Description="yes/no/NA, whether paired local position matches">',
    '##INFO=<ID=MTTRNA_S_PAIR_POS,Number=1,Type=String,Description="Species/source paired genomic position">',
    '##INFO=<ID=MTTRNA_H_PAIR_POS,Number=1,Type=String,Description="Human expected paired genomic position">',
    '##INFO=<ID=MTTRNA_S_PAIR_LIFTED_HPOS,Number=1,Type=String,Description="Species/source paired genomic position lifted to human genomic position through alignment posmap">',
    '##INFO=<ID=MTTRNA_PAIR_POS_MATCH,Number=1,Type=String,Description="yes/no/NA, whether lifted source paired position equals human expected paired position">',
    '##INFO=<ID=MTTRNA_STRICT_MATCH,Number=1,Type=String,Description="yes/no, strict tRNA match: loop requires region+element match; stem requires region+element+pair_state+pair_pos match">',
    '##INFO=<ID=MTTRNA_S_COORD_SPACE,Number=1,Type=String,Description="species tRNA lookup coordinate space: original or rotated">',
    '##INFO=<ID=MTTRNA_S_LOOKUP_CHROM,Number=1,Type=String,Description="species tRNA lookup chrom key (or * when chrom is ignored)">',
    '##INFO=<ID=MTTRNA_S_LOOKUP_POS,Number=1,Type=String,Description="species tRNA lookup position key">',
]


def compare_values(a: str, b: str) -> str:
    if a in {NA, "", "."} or b in {NA, "", "."}:
        return NA
    return "yes" if str(a) == str(b) else "no"


def first_row(rows: List[Dict[str, str]]) -> Optional[Dict[str, str]]:
    return rows[0] if rows else None


def trna_strict_match(info: Dict[str, object]) -> str:
    s_class = str(info.get("MTTRNA_S_CLASS", NA))
    h_class = str(info.get("MTTRNA_H_CLASS", NA))
    region_match = str(info.get("MTTRNA_REGION_MATCH", NA))
    element_match = str(info.get("MTTRNA_ELEMENT_MATCH", NA))
    pair_state_match = str(info.get("MTTRNA_PAIR_STATE_MATCH", NA))
    pair_pos_match = str(info.get("MTTRNA_PAIR_POS_MATCH", NA))
    if s_class == "loop" and h_class == "loop":
        return "yes" if region_match == "yes" and element_match == "yes" else "no"
    if s_class == "stem" and h_class == "stem":
        return "yes" if (region_match == "yes" and element_match == "yes" and pair_state_match == "yes" and pair_pos_match == "yes") else "no"
    return "no"


def annotate_trna_record(
    parts: List[str],
    species_idx: Dict[Tuple[str, int], List[Dict[str, str]]],
    human_idx: Dict[Tuple[str, int], List[Dict[str, str]]],
    species_chrom_norm: str,
    human_chrom_norm: str,
    species_trna_coord_space: str,
    species_lookup_ignore_chrom: bool,
    human_lookup_ignore_chrom: bool,
    map_dict: Dict[int, Tuple[int, str, str, str]],
    p_init: int,
    p_len: int,
    human_len: int,
    human_offset: int,
) -> Dict[str, object]:
    info = parse_info(parts[7])
    sp_chrom = normalize_chrom(info.get("MTLIFT_ORIG_CHROM"), species_chrom_norm)
    coord_space = str(species_trna_coord_space or "original").strip().lower()
    if coord_space not in {"original", "rotated"}:
        coord_space = "original"
    sp_pos = parse_optional_int(info.get("MTLIFT_ORIG_ROT_POS")) if coord_space == "rotated" else parse_optional_int(info.get("MTLIFT_ORIG_POS"))
    hu_chrom = normalize_chrom(parts[0], human_chrom_norm)
    hu_pos = int(parts[1])
    ann: Dict[str, object] = {}
    sp_lookup_chrom = "*" if species_lookup_ignore_chrom else sp_chrom
    ann["MTTRNA_S_COORD_SPACE"] = coord_space
    ann["MTTRNA_S_LOOKUP_CHROM"] = sp_lookup_chrom if sp_lookup_chrom is not None else NA
    ann["MTTRNA_S_LOOKUP_POS"] = sp_pos if sp_pos is not None else NA
    if sp_lookup_chrom is None or sp_pos is None:
        ann["MTTRNA_STATUS"] = "MISSING_SPECIES_COORD"
        return ann
    sp_rows = species_idx.get((sp_lookup_chrom, sp_pos), [])
    hu_lookup_chrom = "*" if human_lookup_ignore_chrom else hu_chrom
    hu_rows = human_idx.get((hu_lookup_chrom, hu_pos), [])
    sp = first_row(sp_rows)
    hu = first_row(hu_rows)
    if sp is None and hu is None:
        ann["MTTRNA_STATUS"] = "NO_SPECIES_OR_HUMAN_TRNA"
    elif sp is None:
        ann["MTTRNA_STATUS"] = "NO_SPECIES_TRNA"
    elif hu is None:
        ann["MTTRNA_STATUS"] = "NO_HUMAN_TRNA"
    else:
        ann["MTTRNA_STATUS"] = "OK"
    def g(row, key):
        return row.get(key, NA) if row else NA
    ann.update({
        "MTTRNA_S_ID": g(sp, "trna_id"),
        "MTTRNA_H_ID": g(hu, "trna_id"),
        "MTTRNA_S_LOCAL": g(sp, "local_pos"),
        "MTTRNA_H_LOCAL": g(hu, "local_pos"),
        "MTTRNA_S_CLASS": g(sp, "struct_class"),
        "MTTRNA_H_CLASS": g(hu, "struct_class"),
        "MTTRNA_REGION_MATCH": compare_values(g(sp, "struct_class"), g(hu, "struct_class")),
        "MTTRNA_S_ELEMENT": g(sp, "struct_element"),
        "MTTRNA_H_ELEMENT": g(hu, "struct_element"),
        "MTTRNA_ELEMENT_MATCH": compare_values(g(sp, "struct_element"), g(hu, "struct_element")),
        "MTTRNA_S_PAIR_TYPE": g(sp, "pair_type"),
        "MTTRNA_H_PAIR_TYPE": g(hu, "pair_type"),
        "MTTRNA_PAIR_TYPE_MATCH": compare_values(g(sp, "pair_type"), g(hu, "pair_type")),
        "MTTRNA_S_PAIR_STATE": g(sp, "pair_state"),
        "MTTRNA_H_PAIR_STATE": g(hu, "pair_state"),
        "MTTRNA_PAIR_STATE_MATCH": compare_values(g(sp, "pair_state"), g(hu, "pair_state")),
        "MTTRNA_S_PAIR_LOCAL": g(sp, "paired_local_pos"),
        "MTTRNA_H_PAIR_LOCAL": g(hu, "paired_local_pos"),
        "MTTRNA_PAIR_LOCAL_MATCH": compare_values(g(sp, "paired_local_pos"), g(hu, "paired_local_pos")),
        "MTTRNA_S_PAIR_POS": g(sp, "paired_genomic_pos"),
        "MTTRNA_H_PAIR_POS": g(hu, "paired_genomic_pos"),
        "MTTRNA_S_PAIR_LIFTED_HPOS": NA,
        "MTTRNA_PAIR_POS_MATCH": NA,
    })
    sp_pair_pos = parse_optional_int(g(sp, "paired_genomic_pos"))
    hu_pair_pos = parse_optional_int(g(hu, "paired_genomic_pos"))
    if sp_pair_pos is not None:
        sp_pair_rot = sp_pair_pos if coord_space == "rotated" else rotate_pos(sp_pair_pos, p_init, p_len)
        hit = map_dict.get(sp_pair_rot)
        if hit is not None:
            h_pair_rot, _usable, _qref, _href = hit
            h_pair_final = restore_human_pos(h_pair_rot, human_offset, human_len)
            ann["MTTRNA_S_PAIR_LIFTED_HPOS"] = h_pair_final
            if hu_pair_pos is not None:
                ann["MTTRNA_PAIR_POS_MATCH"] = "yes" if h_pair_final == hu_pair_pos else "no"
    ann["MTTRNA_STRICT_MATCH"] = trna_strict_match(ann)
    return ann


def annotate_trna_vcf_file(
    in_vcf: str,
    out_vcf: str,
    species_index_path: str,
    human_index_path: str,
    posmap_path: str,
    sample: str,
    sample_to_ref: Dict[str, str],
    ref_rotate_info: Dict[str, Tuple[int, int]],
    human_len: int,
    human_offset: int,
    species_chrom_norm: str,
    human_chrom_norm: str,
    species_trna_coord_space: str,
    species_lookup_ignore_chrom: bool,
    human_lookup_ignore_chrom: bool,
    summary_path: Optional[str] = None,
) -> Counter:
    if sample not in sample_to_ref:
        die(f"sample={sample} not found in sample_ref_file")
    rotate_key = sample_to_ref[sample]
    if rotate_key not in ref_rotate_info:
        die(f"rotate key not found in rotate_pos_file: {rotate_key}")
    p_len, p_init = ref_rotate_info[rotate_key]
    species_idx = load_trna_position_index(species_index_path, chrom_norm=species_chrom_norm, ignore_chrom=species_lookup_ignore_chrom)
    human_idx = load_trna_position_index(human_index_path, chrom_norm=human_chrom_norm, ignore_chrom=human_lookup_ignore_chrom)
    map_dict = build_map_dict(posmap_path)
    stats = Counter()
    with open_text(in_vcf, "rt") as fin, open_text(out_vcf, "wt") as fout:
        header_done = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                for h in TRNA_INFO_LINES:
                    fout.write(h + "\n")
                fout.write(line)
                header_done = True
                break
        if not header_done:
            die(f"Malformed VCF header: {in_vcf}")
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                stats["malformed"] += 1
                continue
            ann = annotate_trna_record(parts, species_idx, human_idx, species_chrom_norm, human_chrom_norm, species_trna_coord_space, species_lookup_ignore_chrom, human_lookup_ignore_chrom, map_dict, p_init, p_len, human_len, human_offset)
            add_info_to_parts(parts, ann)
            status = str(ann.get("MTTRNA_STATUS", "UNKNOWN"))
            stats[f"status_{status}"] += 1
            if ann.get("MTTRNA_REGION_MATCH") == "yes":
                stats["region_match_yes"] += 1
            if ann.get("MTTRNA_PAIR_STATE_MATCH") == "yes":
                stats["pair_state_match_yes"] += 1
            if ann.get("MTTRNA_PAIR_POS_MATCH") == "yes":
                stats["pair_pos_match_yes"] += 1
            stats["written_records"] += 1
            fout.write("\t".join(parts) + "\n")
    if summary_path:
        append_summary(summary_path, stats, "trna")
    return stats


def trnascan_mode_args(mode: str) -> List[str]:
    mode_map = {
        "euk": ["-E"],
        "bact": ["-B"],
        "arch": ["-A"],
        "general": ["-G"],
        "mito_mammal": ["-M", "mammal"],
        "mito_vert": ["-M", "vert"],
        "organellar": ["-O"],
    }
    if mode not in mode_map:
        die(f"Unsupported tRNAscan mode: {mode}")
    return mode_map[mode]


def run_trnascan(fasta: str, prefix: str, mode: str = "mito_mammal", threads: int = 4, trnascan_bin: str = "tRNAscan-SE", extra_args: str = "") -> Tuple[str, str]:
    if shutil.which(trnascan_bin) is None and not Path(trnascan_bin).exists():
        die(f"tRNAscan-SE executable not found: {trnascan_bin}")
    ensure_dir(Path(prefix).parent)
    out = f"{prefix}.trnascan.out"
    ss = f"{prefix}.trnascan.ss"
    stats = f"{prefix}.trnascan.stats"
    bed = f"{prefix}.trnascan.bed"
    fa = f"{prefix}.trnascan.fa"
    cmd = [trnascan_bin] + trnascan_mode_args(mode) + ["-Q", "--thread", str(threads), "-o", out, "-f", ss, "-m", stats, "-b", bed, "-a", fa]
    if extra_args:
        cmd.extend(extra_args.split())
    cmd.append(fasta)
    log("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out, ss


def species_trna_paths(cfg: PipeConfig, sample: str, dirs: Dict[str, Path]) -> Tuple[str, str, str]:
    out_tmpl = cfg.path_get("species_trnascan_out_template", "", sample=sample)
    ss_tmpl = cfg.path_get("species_trnascan_ss_template", "", sample=sample)
    index_tmpl = cfg.path_get("species_trna_index_template", "{trna_index_dir}/{sample}.trna_position_index.tsv.gz", sample=sample)
    if out_tmpl and ss_tmpl:
        return out_tmpl, ss_tmpl, index_tmpl
    prefix = cfg.path_get("species_trnascan_prefix_template", "{outdir}/trnascan/{sample}", sample=sample)
    return f"{prefix}.trnascan.out", f"{prefix}.trnascan.ss", index_tmpl


def ensure_trna_indexes(cfg: PipeConfig, sample: str, dirs: Dict[str, Path]) -> Tuple[str, str]:
    # Human index
    human_index = cfg.path_get("human_trna_index", "{trna_index_dir}/human.trna_position_index.tsv.gz")
    human_out = cfg.path_get("human_trnascan_out", "")
    human_ss = cfg.path_get("human_trnascan_ss", "")
    human_norm = cfg.setting("human_trna_chrom_norm", "none")
    if not human_out or not human_ss:
        prefix = cfg.path_get("human_trnascan_prefix", "{outdir}/trnascan/human")
        human_out = f"{prefix}.trnascan.out"
        human_ss = f"{prefix}.trnascan.ss"
    if not Path(human_index).exists():
        if not Path(human_out).exists() or not Path(human_ss).exists():
            if cfg.setting_bool("run_trnascan_if_missing", False):
                human_fasta = cfg.path_get("human_trna_fasta", cfg.path_get("ref_human_fasta"))
                prefix = cfg.path_get("human_trnascan_prefix", "{outdir}/trnascan/human")
                run_trnascan(human_fasta, prefix, cfg.setting("trnascan_mode", "mito_mammal"), cfg.setting_int("trnascan_threads", 4), cfg.setting("trnascan_bin", "tRNAscan-SE"), cfg.setting("trnascan_extra_args", ""))
                human_out = f"{prefix}.trnascan.out"
                human_ss = f"{prefix}.trnascan.ss"
            else:
                die(f"Human tRNAscan files missing: {human_out}, {human_ss}. Set run_trnascan_if_missing=1 or provide paths.")
        build_trna_position_index("human", human_out, human_ss, human_index, chrom_norm=human_norm)
    # Species index
    sp_out, sp_ss, sp_index = species_trna_paths(cfg, sample, dirs)
    sp_norm = cfg.setting("species_trna_chrom_norm", "none")
    if not Path(sp_index).exists():
        if not Path(sp_out).exists() or not Path(sp_ss).exists():
            if cfg.setting_bool("run_trnascan_if_missing", False):
                fasta = cfg.path_get("species_trna_fasta_template", "", sample=sample)
                if (not fasta) or (not Path(fasta).exists()):
                    fasta = sample_fasta(cfg, sample)
                prefix = cfg.path_get("species_trnascan_prefix_template", "{outdir}/trnascan/{sample}", sample=sample)
                run_trnascan(fasta, prefix, cfg.setting("trnascan_mode", "mito_mammal"), cfg.setting_int("trnascan_threads", 4), cfg.setting("trnascan_bin", "tRNAscan-SE"), cfg.setting("trnascan_extra_args", ""))
                sp_out, sp_ss, sp_index = species_trna_paths(cfg, sample, dirs)
            else:
                die(f"Species tRNAscan files missing for {sample}: {sp_out}, {sp_ss}. Set run_trnascan_if_missing=1 or provide paths.")
        build_trna_position_index(sample, sp_out, sp_ss, sp_index, chrom_norm=sp_norm)
    return sp_index, human_index


def trna_stage(cfg: PipeConfig, sample: str, input_vcfs: List[str], posmap_path: str, summary_path: str, dirs: Dict[str, Path]) -> List[str]:
    species_index, human_index = ensure_trna_indexes(cfg, sample, dirs)
    sample_to_ref = load_sample_ref_map(cfg.path_get("sample_ref_file"))
    ref_rotate_info = load_rotate_info(cfg.path_get("rotate_pos_file"))
    human_len = cfg.setting_int("human_len", 16569)
    human_offset = cfg.setting_int("human_restore_offset", 1325)
    sp_norm = cfg.setting("species_vcf_chrom_norm", cfg.setting("species_trna_chrom_norm", "none"))
    hu_norm = cfg.setting("human_vcf_chrom_norm", cfg.setting("human_trna_chrom_norm", "none"))
    out_dir = ensure_dir(dirs["vcf_trna"] / sample)
    outputs = []
    for in_vcf in input_vcfs:
        stem = strip_vcf_suffix(in_vcf)
        if stem.endswith(".lifted.codon"):
            stem = stem[: -len(".lifted.codon")]
            out_vcf = str(out_dir / f"{stem}.lifted.codon.trna.vcf")
        elif stem.endswith(".lifted.raw"):
            stem = stem[: -len(".lifted.raw")]
            out_vcf = str(out_dir / f"{stem}.lifted.trna.vcf")
        else:
            out_vcf = str(out_dir / f"{stem}.trna.vcf")
        log(f"tRNA annotate: {in_vcf} -> {out_vcf}")
        annotate_trna_vcf_file(in_vcf, out_vcf, species_index, human_index, posmap_path, sample, sample_to_ref, ref_rotate_info, human_len, human_offset, sp_norm, hu_norm, cfg.setting("species_trna_coord_space", "original"), cfg.setting_bool("species_trna_lookup_ignore_chrom", False), cfg.setting_bool("human_trna_lookup_ignore_chrom", False), summary_path)
        outputs.append(out_vcf)
    return outputs


# -----------------------------------------------------------------------------
# Optional final filter
# -----------------------------------------------------------------------------

def record_passes_filter(info: Dict[str, str], mode: str) -> bool:
    mode = mode.strip().lower()
    if mode in {"none", "", "off"}:
        return True
    if mode == "region_policy":
        # 1) coding: require strict codon PASS
        # 2) tRNA: require region/pair-state/pair-position match
        # 3) other noncoding: keep by default
        codon_status = str(info.get("MTCODON_STATUS", NA))
        if codon_status == "PASS":
            return True
        if codon_status not in {NA, "SKIPPED_NONCODING", "MISSING_COORD"}:
            return False

        trna_status = str(info.get("MTTRNA_STATUS", NA))
        is_trna_variant = trna_status in {"OK", "NO_SPECIES_TRNA", "NO_HUMAN_TRNA"}
        if is_trna_variant:
            return info.get("MTTRNA_STRICT_MATCH") == "yes"
        return True
    if mode == "trna_loose_match":
        return (
            info.get("MTTRNA_REGION_MATCH") == "yes" or
            info.get("MTTRNA_PAIR_STATE_MATCH") == "yes" or
            info.get("MTTRNA_PAIR_POS_MATCH") == "yes"
        )
    if mode == "trna_strict_match":
        return info.get("MTTRNA_STRICT_MATCH") == "yes"
    if mode == "codon_pass":
        return info.get("MTCODON_STATUS") == "PASS"
    if mode == "trna_region_match":
        return info.get("MTTRNA_REGION_MATCH") == "yes"
    if mode == "trna_pair_state_match":
        return info.get("MTTRNA_PAIR_STATE_MATCH") == "yes"
    if mode == "trna_pair_pos_match":
        return info.get("MTTRNA_PAIR_POS_MATCH") == "yes"
    if mode == "codon_or_trna":
        return (
            info.get("MTCODON_STATUS") == "PASS" or
            info.get("MTTRNA_REGION_MATCH") == "yes" or
            info.get("MTTRNA_PAIR_STATE_MATCH") == "yes" or
            info.get("MTTRNA_PAIR_POS_MATCH") == "yes"
        )
    die(f"Unsupported final_filter_mode: {mode}")
    return False


def filter_vcf_file(in_vcf: str, out_vcf: str, mode: str, summary_path: Optional[str] = None) -> Counter:
    stats = Counter()
    with open_text(in_vcf, "rt") as fin, open_text(out_vcf, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if not line.strip():
                continue
            stats["input_records"] += 1
            parts = line.rstrip("\n").split("\t")
            info = parse_info(parts[7]) if len(parts) > 7 else {}
            if record_passes_filter(info, mode):
                fout.write(line)
                stats["written_records"] += 1
            else:
                stats["filtered_records"] += 1
    if summary_path:
        append_summary(summary_path, stats, "final_filter")
    return stats


def final_filter_stage(cfg: PipeConfig, sample: str, input_vcfs: List[str], summary_path: str, dirs: Dict[str, Path]) -> List[str]:
    mode = cfg.setting("final_filter_mode", "none")
    out_dir = ensure_dir(dirs["vcf_final"] / sample)
    outputs = []
    for in_vcf in input_vcfs:
        stem = strip_vcf_suffix(in_vcf)
        out_vcf = str(out_dir / f"{stem}.filtered_{mode}.vcf")
        log(f"Final filter ({mode}): {in_vcf} -> {out_vcf}")
        filter_vcf_file(in_vcf, out_vcf, mode, summary_path)
        outputs.append(out_vcf)
    return outputs


# -----------------------------------------------------------------------------
# Sample-level orchestration
# -----------------------------------------------------------------------------

def existing_raw_vcfs(cfg: PipeConfig, sample: str, dirs: Dict[str, Path]) -> List[str]:
    d = dirs["vcf_raw"] / sample
    return sorted(glob.glob(str(d / "*.lifted.raw.vcf")) + glob.glob(str(d / "*.lifted.raw.vcf.gz")))


def run_sample(cfg: PipeConfig, sample: str) -> None:
    dirs = make_outdirs(cfg)
    fasta = sample_fasta(cfg, sample)
    work_dir = str(dirs["tmp"] / sample)
    ensure_dir(work_dir)
    aln_fasta = str(dirs["alignments"] / f"{sample}.aligned.fa")
    if cfg.setting_bool("run_align", True):
        run_alignment(cfg, sample, fasta, aln_fasta, work_dir)
    elif not Path(aln_fasta).exists():
        die(f"run_align=0 but alignment missing: {aln_fasta}")
    if cfg.setting_bool("run_posmap", True) or not Path(dirs["maps"] / f"{sample}.posmap.tsv.gz").exists():
        _cols_path, posmap_path, summary_path = build_posmap_stage(cfg, sample, aln_fasta, dirs)
    else:
        posmap_path = str(dirs["maps"] / f"{sample}.posmap.tsv.gz")
        summary_path = str(dirs["reports"] / f"{sample}.mapping.summary.tsv")
    current_vcfs: List[str] = []
    if cfg.setting_bool("run_vcf_liftover", True) or cfg.setting_bool("run_cov_liftover", False):
        current_vcfs = liftover_stage(cfg, sample, posmap_path, summary_path, dirs)
    else:
        current_vcfs = existing_raw_vcfs(cfg, sample, dirs)
    if cfg.setting_bool("run_codon_annotate", False):
        if not current_vcfs:
            current_vcfs = existing_raw_vcfs(cfg, sample, dirs)
        current_vcfs = codon_stage(cfg, sample, current_vcfs, summary_path, dirs)
    if cfg.setting_bool("run_trna_annotate", False):
        if not current_vcfs:
            current_vcfs = existing_raw_vcfs(cfg, sample, dirs)
        current_vcfs = trna_stage(cfg, sample, current_vcfs, posmap_path, summary_path, dirs)
    if cfg.setting_bool("run_final_filter", False):
        current_vcfs = final_filter_stage(cfg, sample, current_vcfs, summary_path, dirs)
    if cfg.setting_bool("keep_tmp", True) is False:
        shutil.rmtree(work_dir, ignore_errors=True)
    log(f"Done sample={sample}")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Config-driven mtDNA liftover + codon + tRNA annotation pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("list-samples", help="List samples discovered from config")
    p.add_argument("--config", required=True)

    p = sub.add_parser("run-sample", help="Run pipeline for one species/sample")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)

    p = sub.add_parser("run-array", help="Run one array task; uses SLURM_ARRAY_TASK_ID unless --task-id is supplied")
    p.add_argument("--config", required=True)
    p.add_argument("--task-id", type=int, default=None)

    p = sub.add_parser("run-trnascan", help="Run tRNAscan-SE for human or one species according to config")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True, help="Use 'human' for human reference, otherwise species/sample key")

    p = sub.add_parser("build-trna-index", help="Build tRNA position index for human or one species")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True, help="Use 'human' for human reference, otherwise species/sample key")

    p = sub.add_parser("liftover-vcf", help="Run liftover-only for one explicit VCF")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--dropped", required=True)

    p = sub.add_parser("annotate-codon", help="Annotate one lifted VCF with codon match")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)

    p = sub.add_parser("annotate-trna", help="Annotate one lifted VCF with tRNA match")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)

    p = sub.add_parser("filter-vcf", help="Filter one annotated VCF by INFO fields")
    p.add_argument("--config", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--mode", default=None)
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    cfg = PipeConfig(args.config)
    if args.cmd == "list-samples":
        for s in discover_samples(cfg):
            print(s)
        return 0
    if args.cmd == "run-sample":
        run_sample(cfg, args.sample)
        return 0
    if args.cmd == "run-array":
        samples = discover_samples(cfg)
        tid = args.task_id
        if tid is None:
            tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
        if tid < 0 or tid >= len(samples):
            log(f"Task {tid} outside sample count {len(samples)}; skip")
            return 0
        run_sample(cfg, samples[tid])
        return 0
    if args.cmd == "run-trnascan":
        dirs = make_outdirs(cfg)
        if args.sample == "human":
            fasta = cfg.path_get("human_trna_fasta", cfg.path_get("ref_human_fasta"))
            prefix = cfg.path_get("human_trnascan_prefix", "{outdir}/trnascan/human")
        else:
            fasta = cfg.path_get("species_trna_fasta_template", "", sample=args.sample)
            if (not fasta) or (not Path(fasta).exists()):
                fasta = sample_fasta(cfg, args.sample)
            prefix = cfg.path_get("species_trnascan_prefix_template", "{outdir}/trnascan/{sample}", sample=args.sample)
        run_trnascan(fasta, prefix, cfg.setting("trnascan_mode", "mito_mammal"), cfg.setting_int("trnascan_threads", 4), cfg.setting("trnascan_bin", "tRNAscan-SE"), cfg.setting("trnascan_extra_args", ""))
        return 0
    if args.cmd == "build-trna-index":
        dirs = make_outdirs(cfg)
        if args.sample == "human":
            out = cfg.path_get("human_trnascan_out")
            ss = cfg.path_get("human_trnascan_ss")
            index = cfg.path_get("human_trna_index", "{trna_index_dir}/human.trna_position_index.tsv.gz")
            norm = cfg.setting("human_trna_chrom_norm", "none")
            build_trna_position_index("human", out, ss, index, norm)
        else:
            out, ss, index = species_trna_paths(cfg, args.sample, dirs)
            norm = cfg.setting("species_trna_chrom_norm", "none")
            build_trna_position_index(args.sample, out, ss, index, norm)
        return 0
    if args.cmd == "liftover-vcf":
        dirs = make_outdirs(cfg)
        aln_fasta = str(dirs["alignments"] / f"{args.sample}.aligned.fa")
        if not Path(dirs["maps"] / f"{args.sample}.posmap.tsv.gz").exists():
            build_posmap_stage(cfg, args.sample, aln_fasta, dirs)
        posmap_path = str(dirs["maps"] / f"{args.sample}.posmap.tsv.gz")
        map_dict = build_map_dict(posmap_path)
        human_ref_seq = first_fasta_record(cfg.path_get("ref_human_fasta"))[1]
        sample_to_ref = load_sample_ref_map(cfg.path_get("sample_ref_file"))
        ref_rotate_info = load_rotate_info(cfg.path_get("rotate_pos_file"))
        liftover_vcf_file(args.sample, args.sample, args.vcf, args.output, args.dropped, map_dict, sample_to_ref, ref_rotate_info, human_ref_seq, cfg.setting_int("human_len", 16569), cfg.setting_int("human_restore_offset", 1325), cfg.setting("target_chrom", "chrM"), dirs["debug"])
        return 0
    if args.cmd == "annotate-codon":
        primate_lookup = load_all_primate_position_codon_table(cfg.path_get("all_primate_position_codon_table"))
        human_lookup = load_human_codon_lookup(cfg.path_get("human_codon_table"))
        annotate_codon_vcf_file(args.input, args.output, args.sample, primate_lookup, human_lookup, cfg.setting_bool("strict_phase_match", True))
        return 0
    if args.cmd == "annotate-trna":
        dirs = make_outdirs(cfg)
        sp_idx, hu_idx = ensure_trna_indexes(cfg, args.sample, dirs)
        annotate_trna_vcf_file(args.input, args.output, sp_idx, hu_idx, str(dirs["maps"] / f"{args.sample}.posmap.tsv.gz"), args.sample, load_sample_ref_map(cfg.path_get("sample_ref_file")), load_rotate_info(cfg.path_get("rotate_pos_file")), cfg.setting_int("human_len", 16569), cfg.setting_int("human_restore_offset", 1325), cfg.setting("species_vcf_chrom_norm", cfg.setting("species_trna_chrom_norm", "none")), cfg.setting("human_vcf_chrom_norm", cfg.setting("human_trna_chrom_norm", "none")), cfg.setting("species_trna_coord_space", "original"), cfg.setting_bool("species_trna_lookup_ignore_chrom", False), cfg.setting_bool("human_trna_lookup_ignore_chrom", False))
        return 0
    if args.cmd == "filter-vcf":
        mode = args.mode if args.mode is not None else cfg.setting("final_filter_mode", "none")
        filter_vcf_file(args.input, args.output, mode)
        return 0
    die(f"Unknown command: {args.cmd}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
