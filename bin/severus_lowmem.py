#!/usr/bin/env python3
"""
severus_lowmem.py - memory-bounded drop-in replacement for the `severus` command.

    severus_lowmem.py [Severus options ...] [lowmem options ...]

All Severus options are accepted unchanged and all Severus outputs are written as
usual. Extra options (run with --lowmem-help to list them):

    --tmp-dir PATH            directory for scratch files; a private subdirectory is
                              created there and removed at the end [<out-dir>]
    --no-prescreen            hand every read to the breakpoint caller (stock behaviour)
    --prescreen-window INT    breakend bin for junction recurrence, bp [500]
    --prescreen-min-reads INT reads (tumour+normal pooled) that must share a junction [2]
    --buckets INT             read-name hash buckets per BAM [auto, ~2M segments each]
    --keep-tmp                keep the scratch directory

What it changes
---------------
Stock Severus parses every alignment of every BAM into Python objects and keeps all of
them in memory while it computes read statistics, the background mismatch rate,
mismatch-rate windows and coverage histograms, and then passes all of them to the
breakpoint caller. Memory therefore scales with the number of aligned segments.

This driver runs the installed Severus (severus.main.main) and replaces only the read
loading / annotation stage (get_all_reads_parallel, update_segments_by_read and the
segment part of update_coverage_hist):

1. Each 10 Mb chunk (the same chunks Severus uses) is parsed in a worker with Severus's
   own get_segment(). Unsplit reads become Severus's compact 9-column rows; segments are
   written to disk, hash-partitioned by read name.
2. Reads are regrouped bucket by bucket and passed through Severus's own
   remove_dedup_segments, resolve_read_vntr and label_reads. Read N50/N90, the
   background mismatch rate, mismatch-rate windows, coverage histograms and read_qual
   counts are accumulated exactly as Severus computes them (including its chunking).
3. Severus's breakpoint caller and graph builder then run unchanged on the reads that
   can contribute to calls.

Prescreen (on by default)
-------------------------
A read whose segments come from two or more separate alignments (primary +
supplementary) is passed to the breakpoint caller only if at least one of its junctions
(any pair of its segments, after VNTR resolution) is shared, within +/- one bin at both
breakends, by another read in any of the input BAMs. Other split reads still contribute
to all statistics and coverage, and their insertion segments and clipped ends are still
passed on; only their junctions are withheld. Reads split only within one alignment
(CIGAR deletions/insertions) and unsplit reads are always passed on.

A junction seen in a single read can never reach Severus's minimum support, but in
stock Severus such reads still join the 50 bp breakend clusters, where they can shift
breakpoint positions or mark real breakpoints as single/FAIL_CONN_CONS. With chimeric
libraries they also dominate memory. With --no-prescreen the output is the same as
stock Severus; memory is then bounded for parsing but not for breakpoint calling.

Resources
---------
Memory is about 1 GB per thread (one 10 Mb chunk or one read bucket per worker) plus
the reads passed to the caller, independent of BAM size. Scratch space is about 80
bytes per aligned segment (~10% of the BAM size) and is deleted at the end. Progress,
worker memory and scratch use are logged to severus.log as "[lowmem]" lines.

Tested with Severus 1.6 and 1.7. Requires pysam and numpy (both Severus dependencies).
"""

__version__ = "0.1.0"

import argparse
import hashlib
import inspect
import logging
import multiprocessing
import os
import pickle
import resource
import shutil
import sys
import tempfile
import time
import zlib
from collections import defaultdict

import numpy as np
import pysam

import severus.main as sv_main
import severus.bam_processing as sv_bp
import severus.resolve_vntr as sv_vntr
import severus.breakpoint_finder as sv_bf
from severus.__version__ import __version__ as SEVERUS_VERSION

logger = logging.getLogger()

TESTED_VERSIONS = ("1.6", "1.7")
CHUNK_SIZE = 10000000              # bam_processing.get_all_reads_parallel
COV_WINDOW = sv_bp.COV_WINDOW
COV_WINDOW_MM = sv_bp.COV_WINDOW_MM
NUM_HAPLOTYPES = sv_bp.NUM_HAPLOTYPES
MM_BITS = 21
MM_OFF = 1 << 20
MM_MASK = (1 << MM_BITS) - 1
REGION_W = 10000                   # mismatch windows per reduction region
GROUP_BATCH_ROWS = 200000          # rows converted to Python objects at a time
TARGET_BUCKET_ROWS = 1500000

FL_INS, FL_CLIP, FL_PRIMARY, FL_ID = 1, 2, 4, 8

SEG_DTYPE = np.dtype([
    ("h1", "<u8"), ("h2", "<u8"), ("cg", "<u4"), ("rec", "<u4"), ("sidx", "<u4"),
    ("tid", "<i4"), ("flags", "<u1"), ("strand", "<i1"), ("mapq", "<i2"),
    ("align_start", "<i8"), ("read_start", "<i8"), ("read_end", "<i8"),
    ("ref_start", "<i8"), ("ref_end", "<i8"), ("ref_start_ori", "<i8"),
    ("ref_end_ori", "<i8"), ("read_length", "<i8"), ("align_len", "<i8"),
    ("segment_length", "<i8"), ("haplotype", "<i8"), ("PS", "<i8"), ("mm", "<i8"),
    ("err", "<i8"), ("rec_s", "<i8"), ("rec_e", "<i8"), ("obj", "<i4")])

INT_FIELDS = ("align_start", "read_start", "read_end", "ref_start", "ref_end",
              "ref_start_ori", "ref_end_ori", "read_length", "align_len",
              "segment_length", "haplotype", "PS")
SLOT_FROM_COL = {"mismatch_rate": "mm", "error_rate": "err"}
HANDLED_SLOTS = set(INT_FIELDS) | {
    "mismatch_rate", "error_rate", "mapq", "strand", "read_id", "ref_id", "genome_id",
    "is_pass", "is_insertion", "is_clipped", "is_primary", "ins_seq", "ins_pos",
    "bp_pos", "ID"}

RS = sv_bp.ReadSegment
RS_SLOTS = tuple(RS.__slots__)
HAS_PS = "PS" in RS_SLOTS
INT_FIELDS_USED = tuple(f for f in INT_FIELDS if f in RS_SLOTS)
EXTRA_SLOTS = tuple(s for s in RS_SLOTS if s not in HANDLED_SLOTS)
_BIG = 1 << 62

CTX = {}       # driver state in the parent process
_G = {}        # read-only parameters inherited by forked workers


class ResidualRead(list):
    """Insertion and clipped-end segments of a split read whose junctions were withheld
    by the prescreen. Clipped ends were already extracted, so extract_clipped_end skips it."""
    __slots__ = ()


# ----------------------------------------------------------------------------- utils

def _fork_pool(n):
    return multiprocessing.get_context("fork").Pool(n)


def _rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576


def _dir_gb(path):
    tot = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                tot += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return tot / 1e9


def _name_hash(name):
    d = hashlib.blake2b(name.encode(), digest_size=16).digest()
    return int.from_bytes(d[:8], "little"), int.from_bytes(d[8:], "little")


def _merge_hist(acc, vals, cnts):
    for v, c in zip(vals.tolist(), cnts.tolist()):
        acc[v] = acc.get(v, 0) + c


def _unique_counts(a):
    if len(a) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    v, c = np.unique(np.asarray(a, dtype=np.int64), return_counts=True)
    return v, c.astype(np.int64)


def _weighted_unique(vals, w):
    keep = w > 0
    vals, w = vals[keep], w[keep]
    if len(vals) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    u, inv = np.unique(vals, return_inverse=True)
    s = np.zeros(len(u), np.int64)
    np.add.at(s, inv, w)
    return u, s


def _hist_arrays(d):
    if not d:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    k = np.array(sorted(d), dtype=np.int64)
    return k, np.array([d[x] for x in k.tolist()], dtype=np.int64)


def _quantile_from_hist(vals, cnts, q):
    """np.quantile(x, q) (linear method) for the multiset given as sorted values+counts."""
    n = int(cnts.sum())
    virtual = np.float64(n - 1) * np.float64(q)
    prev = np.floor(virtual)
    gamma = float(virtual - prev)
    i0 = min(max(int(prev), 0), n - 1)
    i1 = min(i0 + 1, n - 1)
    cum = np.cumsum(cnts)
    v0 = vals[np.searchsorted(cum, i0, side="right")]
    v1 = vals[np.searchsorted(cum, i1, side="right")]
    return np.quantile(np.array([v0, v1], dtype=np.int64), gamma)


def _calc_nx_from_hist(vals, cnts, norm_len, rate):
    """bam_processing._calc_nx on a histogram (vals ascending)."""
    target = rate * np.int64(norm_len)
    s = 0
    for v, c in zip(vals[::-1].tolist(), cnts[::-1].tolist()):
        if c and s + v * c > target:
            return np.int64(v)
        s += v * c
    return 0


def _expand_ranges(starts, ends, vals):
    """All (i, val) for i in [start, end] inclusive, per range."""
    lens = ends - starts + 1
    keep = lens > 0
    starts, lens, vals = starts[keep], lens[keep], vals[keep]
    if len(starts) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    tot = int(lens.sum())
    rep_start = np.repeat(starts, lens)
    first = np.repeat(np.cumsum(lens) - lens, lens)
    idx = rep_start + (np.arange(tot, dtype=np.int64) - first)
    return idx, np.repeat(vals, lens)


def _mm_keys(gw, vals):
    vals = np.clip(vals, -MM_OFF, MM_OFF - 1)
    return (gw.astype(np.int64) << MM_BITS) | (vals.astype(np.int64) + MM_OFF)


def _save_region_sorted(path, keys, n_regions):
    keys = np.sort(keys)
    np.save(path, keys)
    bounds = (np.arange(n_regions + 1, dtype=np.int64) * REGION_W) << MM_BITS
    return np.searchsorted(keys, bounds).astype(np.int64)


# ------------------------------------------------------------ phase 1: parse chunks

def _p1_worker(task):
    bam_idx, cg, j, ctg, start, end = task
    P = _G
    genome_id = P["genome_id"]
    sv_size, use_supp, ignore_hp = P["sv_size"], P["use_supp"], P["ignore_hp"]
    nb = P["nbuckets"]
    get_segment = sv_bp.get_segment
    six = P["gs_six"]
    tmp = P["tmp"]

    cap = 1 << 15
    buf = np.zeros(cap, dtype=SEG_DTYPE)
    n = 0
    objs = []
    names = {}
    inf = []
    rec = -1
    seg_mm_vals, seg_mm_w = [], []

    with pysam.AlignmentFile(P["bam_path"], "rb") as bam:
        for aln in bam.fetch(ctg, start, end, multiple_iterators=True):
            if aln.is_secondary or aln.is_unmapped:
                continue
            rec += 1
            if six:
                segs, rinf = get_segment(aln, genome_id, sv_size, use_supp, j, ignore_hp)
            else:
                segs, rinf = get_segment(aln, genome_id, sv_size, use_supp, j)
            if segs:
                qn = aln.query_name
                h1, h2 = _name_hash(qn)
                rname = aln.reference_name
                rs, re_ = aln.reference_start, aln.reference_end
                exp_id = f"{j},{rs},{re_}"
                if not aln.is_supplementary and (h1, h2) not in names:
                    names[(h1, h2)] = qn
                for k, s in enumerate(segs):
                    if n == cap:
                        cap *= 2
                        nbuf = np.zeros(cap, dtype=SEG_DTYPE)
                        nbuf[:n] = buf[:n]
                        buf = nbuf
                    over = {}
                    vals = []
                    for f in INT_FIELDS_USED:
                        v = getattr(s, f)
                        if type(v) is not int or not -_BIG < v < _BIG:
                            over[f] = v
                            v = 0
                        vals.append(v)
                    if not HAS_PS:
                        vals.append(0)
                    mm, err, mq, st = s.mismatch_rate, s.error_rate, s.mapq, s.strand
                    if type(mm) is not int or not -_BIG < mm < _BIG:
                        over["mismatch_rate"] = mm
                        mm = 0
                    if type(err) is not int or not -_BIG < err < _BIG:
                        over["error_rate"] = err
                        err = 0
                    if type(mq) is not int or not 0 <= mq < 32768:
                        over["mapq"] = mq
                        mq = 0
                    if st not in (1, -1) or type(st) is not int:
                        over["strand"] = st
                        st = 1
                    fl = 0
                    if s.is_insertion is True:
                        fl |= FL_INS
                    elif s.is_insertion is not False:
                        over["is_insertion"] = s.is_insertion
                    if s.is_clipped is True:
                        fl |= FL_CLIP
                    elif s.is_clipped is not False:
                        over["is_clipped"] = s.is_clipped
                    if s.is_primary is True:
                        fl |= FL_PRIMARY
                    elif s.is_primary is not None:
                        over["is_primary"] = s.is_primary
                    if s.ID == exp_id:
                        fl |= FL_ID
                    elif s.ID is not None:
                        over["ID"] = s.ID
                    if s.ins_seq is not None:
                        over["ins_seq"] = s.ins_seq
                    if s.ins_pos is not None:
                        over["ins_pos"] = s.ins_pos
                    if s.bp_pos is not None:
                        over["bp_pos"] = s.bp_pos
                    if s.is_pass != "":
                        over["is_pass"] = s.is_pass
                    if s.read_id != qn:
                        over["read_id"] = s.read_id
                    if s.ref_id != rname:
                        over["ref_id"] = s.ref_id
                    if s.genome_id != genome_id:
                        over["genome_id"] = s.genome_id
                    for slot in EXTRA_SLOTS:
                        over[slot] = getattr(s, slot, None)
                    oi = -1
                    if over:
                        oi = len(objs)
                        objs.append(over)
                    buf[n] = (h1, h2, cg, rec, k, j, fl, st, mq, *vals, mm, err, rs, re_, oi)
                    n += 1
                    if not s.is_clipped and not s.is_insertion:
                        seg_mm_vals.append(s.mismatch_rate)
                        seg_mm_w.append(max((s.ref_end - s.ref_start) // 2000 + 1, 0))
            if len(rinf):
                inf.append(rinf)

    # read_info rows exactly as bam_processing.get_all_reads returns them
    ni = len(inf)
    if ni == 0:
        rows = np.zeros((0, 9), dtype=np.int64)
    else:
        rows = np.vstack(inf).astype(np.int64)
        if ni % 10000 == 0:
            rows = np.vstack([rows, np.zeros((9999, 9), dtype=np.int64)])
        else:
            rows = rows[:ni - 1]
    del inf

    out = {"bam_idx": bam_idx, "cg": cg, "n_rows": len(rows), "n_seg": n}

    # background mismatch rate contributions (background_mm_rat)
    bv = [rows[:, 6]] if len(rows) else []
    bw = [rows[:, 3] // 2000 + 1] if len(rows) else []
    if seg_mm_vals:
        bv.append(np.asarray(seg_mm_vals, dtype=np.int64))
        bw.append(np.asarray(seg_mm_w, dtype=np.int64))
    if bv:
        out["bgmm"] = _weighted_unique(np.concatenate(bv), np.concatenate(bw))
    else:
        out["bgmm"] = _unique_counts([])

    # read statistics contributions (get_read_statistics)
    out["st_rl"] = _unique_counts(rows[:, 3])
    out["st_al"] = _unique_counts(rows[:, 4])
    out["st_err"] = _unique_counts(rows[:, 7])
    out["st_mm"] = _unique_counts(rows[:, 6])

    # mismatch windows from compact rows (update_mm_hist)
    pref = os.path.join(tmp, f"c{cg:06d}")
    if len(rows):
        c0 = int(rows[0][0])
        ws = rows[:, 1] // COV_WINDOW_MM
        we = rows[:, 2] // COV_WINDOW_MM
        w, v = _expand_ranges(ws, we, rows[:, 6])
        nwin = P["mm_nwin"][c0]
        if (w < -nwin).any() or (w >= nwin).any():
            raise IndexError("mismatch window out of range")
        w = np.where(w < 0, w + nwin, w)
        out["mm_offs"] = _save_region_sorted(pref + ".rmm.npy",
                                             _mm_keys(P["mm_off"][c0] + w, v),
                                             P["n_regions"])
        out["mm_path"] = pref + ".rmm.npy"
        # compact rows for coverage (update_cov_hist): start, end, rl, al, hp, mm, mapq
        np.save(pref + ".rows.npy", rows[:, [1, 2, 3, 4, 5, 6, 8]])
        out["rows_path"] = pref + ".rows.npy"
        out["rows_chr"] = c0
    else:
        out["mm_path"] = None
        out["rows_path"] = None

    # segments, side objects and names, partitioned by read-name hash
    buf = buf[:n]
    bkt = (buf["h1"] % np.uint64(nb)).astype(np.int64)
    order = np.argsort(bkt, kind="stable")
    buf, bkt = buf[order], bkt[order]
    bounds = np.searchsorted(bkt, np.arange(nb + 1))
    name_b = defaultdict(lambda: ([], [], []))
    for (h1, h2), qn in names.items():
        t = name_b[h1 % nb]
        t[0].append(h1)
        t[1].append(h2)
        t[2].append(qn)
    seg_off = np.zeros(nb + 1, np.int64)
    obj_off = np.zeros(nb + 1, np.int64)
    nam_off = np.zeros(nb + 1, np.int64)
    pos = 0
    with open(pref + ".seg", "wb") as fh:
        blobs_obj = []
        for k in range(nb):
            part = buf[bounds[k]:bounds[k + 1]]
            if len(part):
                oi = part["obj"]
                if (oi >= 0).any():
                    part = part.copy()
                    sel = np.flatnonzero(part["obj"] >= 0)
                    blobs_obj.append([objs[i] for i in part["obj"][sel].tolist()])
                    part["obj"][sel] = np.arange(len(sel), dtype=np.int32)
                else:
                    blobs_obj.append(None)
                blob = zlib.compress(part.tobytes(), 1)
                fh.write(blob)
                pos += len(blob)
            else:
                blobs_obj.append(None)
            seg_off[k + 1] = pos
        obj_off[0] = pos
        for k in range(nb):
            if blobs_obj[k] is not None:
                blob = zlib.compress(pickle.dumps(blobs_obj[k], protocol=4), 1)
                fh.write(blob)
                pos += len(blob)
            obj_off[k + 1] = pos
        nam_off[0] = pos
        for k in range(nb):
            if k in name_b:
                h1s, h2s, qns = name_b[k]
                blob = zlib.compress(pickle.dumps(
                    (np.array(h1s, dtype=np.uint64), np.array(h2s, dtype=np.uint64), qns),
                    protocol=4), 1)
                fh.write(blob)
                pos += len(blob)
            nam_off[k + 1] = pos
    out.update(seg_path=pref + ".seg", seg_off=seg_off, obj_off=obj_off, nam_off=nam_off,
               rss=_rss_gb())
    return out


# --------------------------------------------------- bucket processing (per read)

def _load_bucket(k, with_names):
    parts, objs = [], []
    nmap = {}
    for ch in _G["chunk_files"]:
        a, b = int(ch["seg_off"][k]), int(ch["seg_off"][k + 1])
        if b <= a:
            continue
        with open(ch["seg_path"], "rb") as fh:
            fh.seek(a)
            arr = np.frombuffer(zlib.decompress(fh.read(b - a)), dtype=SEG_DTYPE).copy()
            oa, ob = int(ch["obj_off"][k]), int(ch["obj_off"][k + 1])
            if ob > oa:
                fh.seek(oa)
                lst = pickle.loads(zlib.decompress(fh.read(ob - oa)))
                sel = arr["obj"] >= 0
                arr["obj"][sel] += len(objs)
                objs.extend(lst)
            if with_names:
                na, nb_ = int(ch["nam_off"][k]), int(ch["nam_off"][k + 1])
                if nb_ > na:
                    fh.seek(na)
                    h1s, h2s, qns = pickle.loads(zlib.decompress(fh.read(nb_ - na)))
                    for x, y, q in zip(h1s.tolist(), h2s.tolist(), qns):
                        nmap[(x, y)] = q
        parts.append(arr)
    if parts:
        arr = np.concatenate(parts)
    else:
        arr = np.zeros(0, dtype=SEG_DTYPE)
    order = np.lexsort((arr["sidx"], arr["rec"], arr["cg"], arr["h2"], arr["h1"]))
    return arr[order], objs, nmap


def _iter_reads(arr, objs, nmap):
    """Yield (key, n_records, read) with ReadSegment objects rebuilt in Severus order."""
    if len(arr) == 0:
        return
    newg = np.ones(len(arr), dtype=bool)
    newg[1:] = (arr["h1"][1:] != arr["h1"][:-1]) | (arr["h2"][1:] != arr["h2"][:-1])
    gstart = np.flatnonzero(newg)
    gend = np.append(gstart[1:], len(arr))
    refnames = _G["refnames"]
    genome_id = _G["genome_id"]
    new = RS.__new__
    b0 = 0
    while b0 < len(gstart):
        b1 = b0
        rows = 0
        while b1 < len(gstart) and (rows < GROUP_BATCH_ROWS or b1 == b0):
            rows += gend[b1] - gstart[b1]
            b1 += 1
        lo, hi = int(gstart[b0]), int(gend[b1 - 1])
        sub = arr[lo:hi]
        c = {f: sub[f].tolist() for f in SEG_DTYPE.names}
        for g in range(b0, b1):
            s0, s1 = int(gstart[g]) - lo, int(gend[g]) - lo
            h = (c["h1"][s0], c["h2"][s0])
            name = nmap.get(h) if nmap is not None else None
            if name is None:
                name = f"{h[0]:016x}{h[1]:016x}"
            read = []
            recs = set()
            for i in range(s0, s1):
                s = new(RS)
                fl = c["flags"][i]
                tid = c["tid"][i]
                s.align_start = c["align_start"][i]
                s.read_start = c["read_start"][i]
                s.read_end = c["read_end"][i]
                s.ref_start = c["ref_start"][i]
                s.ref_end = c["ref_end"][i]
                s.ref_start_ori = c["ref_start_ori"][i]
                s.ref_end_ori = c["ref_end_ori"][i]
                s.read_id = name
                s.ref_id = refnames[tid]
                s.strand = c["strand"][i]
                s.read_length = c["read_length"][i]
                s.align_len = c["align_len"][i]
                s.segment_length = c["segment_length"][i]
                s.haplotype = c["haplotype"][i]
                s.mapq = c["mapq"][i]
                s.genome_id = genome_id
                s.mismatch_rate = c["mm"][i]
                s.is_pass = ""
                s.is_insertion = bool(fl & FL_INS)
                s.is_clipped = bool(fl & FL_CLIP)
                s.error_rate = c["err"][i]
                s.ins_seq = None
                s.is_primary = True if fl & FL_PRIMARY else None
                s.ins_pos = None
                s.bp_pos = None
                s.ID = f"{tid},{c['rec_s'][i]},{c['rec_e'][i]}" if fl & FL_ID else None
                if HAS_PS:
                    s.PS = c["PS"][i]
                oi = c["obj"][i]
                if oi >= 0:
                    for kk, vv in objs[oi].items():
                        setattr(s, kk, vv)
                read.append(s)
                recs.add((tid, c["rec_s"][i], c["rec_e"][i], c["strand"][i], fl & FL_PRIMARY))
            yield (c["cg"][s0], c["rec"][s0], c["sidx"][s0]), len(recs), read
        b0 = b1


def _dedup_and_vntr(read):
    lst = [read]
    sv_vntr.remove_dedup_segments(lst)
    read = lst[0]
    vntr = _G.get("vntr_list")
    if vntr is not None and read:
        new = sv_vntr.resolve_read_vntr(read, vntr, _G["sv_size"])
        read = new if new else []
    return read


def _junction_codes(read):
    """Junction codes for every ordered pair of non-insertion, non-clipped segments."""
    split = [s for s in read if not s.is_insertion and not s.is_clipped]
    if len(split) < 2:
        return set()
    split.sort(key=lambda x: (x.align_start, x.read_start))
    W, sh, eb = _G["pw"], _G["bin_bits"] + 1, _G["ebits"]
    tidmap = _G["tidmap"]
    ex, en = [], []
    for s in split:
        t = tidmap[s.ref_id] << sh
        if s.strand == 1:
            ex.append(t | ((max(s.ref_end, 0) // W) << 1) | 1)
            en.append(t | ((max(s.ref_start, 0) // W) << 1))
        else:
            ex.append(t | ((max(s.ref_start, 0) // W) << 1))
            en.append(t | ((max(s.ref_end, 0) // W) << 1) | 1)
    out = set()
    m = len(split)
    for a in range(m - 1):
        e1 = ex[a]
        for b in range(a + 1, m):
            e2 = en[b]
            out.add(((e1 << eb) | e2) if e1 <= e2 else ((e2 << eb) | e1))
    return out


def _neighbours(code):
    eb = _G["ebits"]
    lo, hi = code >> eb, code & ((1 << eb) - 1)
    res = set()
    for da in (-2, 0, 2):
        for db in (-2, 0, 2):
            a, b = lo + da, hi + db
            if a < 0 or b < 0:
                continue
            res.add(((a << eb) | b) if a <= b else ((b << eb) | a))
    return res


def _supported(codes):
    U, N = _G["rc_codes"], _G["rc_n"]
    if not codes or len(U) == 0:
        return False
    need = _G["pmin"] - 1
    arr = np.fromiter(codes, dtype=np.int64, count=len(codes))
    idx = np.searchsorted(U, arr)
    idx[idx >= len(U)] = 0
    hit = U[idx] == arr
    if not hit.any():
        return False
    for c, nn in zip(arr[hit].tolist(), N[idx[hit]].tolist()):
        own = len(_neighbours(c) & codes)
        if nn - own >= need:
            return True
    return False


def _b1_worker(k):
    """Bucket pass 1: raw read statistics, mismatch windows from segments, junction codes."""
    arr, objs, _ = _load_bucket(k, False)
    rl0, err0, mm0, alnl = [], [], [], []
    ws, we, wv = [], [], []
    codes = []
    n_reads = 0
    mm_off, mm_nwin, tidmap, reflen = _G["mm_off"], _G["mm_nwin"], _G["tidmap"], _G["reflen"]
    prescreen = _G["prescreen"]
    for key, nrec, read in _iter_reads(arr, objs, None):
        n_reads += 1
        r0 = read[0]
        rl0.append(r0.read_length)
        err0.append(r0.error_rate)
        mm0.append(r0.mismatch_rate)
        alnl.append(sum(s.segment_length for s in read if not s.is_clipped and not s.is_insertion))
        read = _dedup_and_vntr(read)
        for s in read:
            if s.is_clipped or s.is_insertion:
                continue
            t = tidmap[s.ref_id]
            a = s.ref_start_ori // COV_WINDOW_MM
            b = min(s.ref_end_ori, reflen[t]) // COV_WINDOW_MM
            if a < 0 or b < 0:
                nwin = mm_nwin[t]
                for i in range(a, b + 1):
                    if i < -nwin:
                        raise IndexError("mismatch window out of range")
                    ws.append(mm_off[t] + (i if i >= 0 else i + nwin))
                    we.append(ws[-1])
                    wv.append(s.mismatch_rate)
            else:
                ws.append(mm_off[t] + a)
                we.append(mm_off[t] + b)
                wv.append(s.mismatch_rate)
        if prescreen and nrec >= 2:
            codes.extend(_junction_codes(read))
    w, v = _expand_ranges(np.asarray(ws, np.int64), np.asarray(we, np.int64),
                          np.asarray(wv, np.int64))
    path = os.path.join(_G["tmp"], f"b{_G['bam_idx']}_{k:05d}.smm.npy")
    offs = _save_region_sorted(path, _mm_keys(w, v), _G["n_regions"])
    del w, v, ws, we, wv
    codes = np.sort(np.asarray(codes, dtype=np.int64))
    cpath = None
    if len(codes):
        cpath = os.path.join(_G["tmp"], f"b{_G['bam_idx']}_{k:05d}.codes.npy")
        np.save(cpath, codes)
    return {"k": k, "n_reads": n_reads, "n_rows": len(arr),
            "rl": _unique_counts(rl0), "err": _unique_counts(err0), "mm": _unique_counts(mm0),
            "al": _unique_counts(alnl), "codes_path": cpath, "n_codes": len(codes),
            "codes_sample": codes[::997] >> _G["ebits"],
            "mm_path": path, "mm_offs": offs, "rss": _rss_gb()}


def _b2_worker(k):
    """Bucket pass 2: labels, read_qual, coverage, and the reads handed to the caller."""
    arr, objs, nmap = _load_bucket(k, True)
    P = _G
    min_mapq, bg_mm, mm_high = P["min_mapq"], P["bg_mm"], P["mm_hist_high"]
    min_al, multisample, min_len = P["min_aligned_length"], P["multisample"], P["min_len"]
    lr7 = P["lr_seven"]
    label_reads = sv_bp.label_reads
    cov_off, reflen, tidmap = P["cov_off"], P["reflen"], P["tidmap"]
    prescreen = P["prescreen"]
    rq_n, rq_len, rq_first = defaultdict(int), defaultdict(int), {}
    cs = {h: [] for h in range(NUM_HAPLOTYPES)}
    ce = {h: [] for h in range(NUM_HAPLOTYPES)}
    out = []
    n_full = n_res = n_withheld = 0
    for key, nrec, read in _iter_reads(arr, objs, nmap):
        read = _dedup_and_vntr(read)
        if not read:
            continue
        if lr7:
            label_reads(read, min_mapq, bg_mm, mm_high, min_al, multisample, min_len)
        else:
            label_reads(read, min_mapq, bg_mm, mm_high, min_al, multisample)
        for pos, s in enumerate(read):
            if s.is_clipped or s.is_insertion or "vntr_only" in s.is_pass:
                continue
            lab = s.is_pass
            rq_n[lab] += 1
            rq_len[lab] += s.ref_end - s.ref_start
            if lab not in rq_first or (key, pos) < rq_first[lab]:
                rq_first[lab] = (key, pos)
        for s in read:
            if s.is_pass == "PASS" and not s.is_insertion and not s.is_clipped:
                if s.haplotype not in cs:
                    raise KeyError((s.genome_id, s.haplotype, s.ref_id))
                t = tidmap[s.ref_id]
                a = s.ref_start_ori // COV_WINDOW + 1
                b = min(s.ref_end_ori, reflen[t]) // COV_WINDOW
                if b <= a:
                    continue
                if a < 0:
                    ncw = reflen[t] // COV_WINDOW + 1
                    for i in range(a, b):
                        if i < -ncw:
                            raise IndexError("coverage window out of range")
                        g = cov_off[t] + (i if i >= 0 else i + ncw)
                        cs[s.haplotype].append(g)
                        ce[s.haplotype].append(g + 1)
                else:
                    cs[s.haplotype].append(cov_off[t] + a)
                    ce[s.haplotype].append(cov_off[t] + b)
        if prescreen and nrec >= 2:
            split = sum(1 for s in read if not s.is_insertion and not s.is_clipped)
            if split >= 2 and not _supported(_junction_codes(read)):
                n_withheld += 1
                sv_bp.extract_clipped_end([read])
                res = ResidualRead(s for s in read if s.is_insertion or s.is_clipped)
                if res:
                    out.append((key, res))
                    n_res += 1
                continue
        out.append((key, read))
        n_full += 1
    ncov = P["n_cov"] + 1
    cov = {}
    for h in range(NUM_HAPLOTYPES):
        if cs[h]:
            d = np.bincount(np.asarray(cs[h], np.int64), minlength=ncov)[:ncov].astype(np.int32)
            d -= np.bincount(np.asarray(ce[h], np.int64), minlength=ncov)[:ncov].astype(np.int32)
            cov[h] = d
    path = os.path.join(P["tmp"], f"r{P['bam_idx']}_{k:05d}.pkl")
    with open(path, "wb") as fh:
        pickle.dump(out, fh, protocol=4)
    return {"k": k, "rq_n": dict(rq_n), "rq_len": dict(rq_len), "rq_first": rq_first,
            "cov": cov, "path": path, "n_full": n_full, "n_res": n_res,
            "n_withheld": n_withheld, "rss": _rss_gb()}


# ------------------------------------------------ compact-row coverage, mismatch flags

def _ric_worker(task):
    path, c0, cg = task
    rows = np.load(path)
    P = _G
    st, en, rl, al, hp, mm, mq = (rows[:, i] for i in range(7))
    min_mapq, bg, n90 = P["min_mapq"], P["bg_mm_bam"], P["n90_bam"]
    passm = (mq > min_mapq) & (mm < bg) & (rl > n90) & (np.abs(rl - al) < al)
    rest = ~passm
    f_mapq = rest & (mq < min_mapq)
    rest &= ~f_mapq
    f_mm = rest & (mm > bg)
    rest &= ~f_mm
    f_rl = rest & (rl < n90)
    rest &= ~f_rl
    f_al = rest & (np.abs(rl - al) > al)
    labs = {}
    for lab, m in (("PASS", passm), ("FAIL_MAPQ", f_mapq), ("FAIL_MM", f_mm),
                   ("FAIL_READLEN", f_rl), ("FAIL_ALNLEN", f_al)):
        if m.any():
            labs[lab] = (int(m.sum()), int(al[m].sum()), (cg, int(np.argmax(m))))
    cov = {}
    if passm.any():
        bad = set(np.unique(hp[passm]).tolist()) - set(range(NUM_HAPLOTYPES))
        if bad:
            raise KeyError(("coverage haplotype", sorted(bad)))
        ws = st[passm] // COV_WINDOW
        we = en[passm] // COV_WINDOW
        w0 = int(ws.min())
        L = int(we.max()) - w0 + 2
        for h in range(NUM_HAPLOTYPES):
            m = hp[passm] == h
            if m.any():
                d = np.bincount(ws[m] - w0, minlength=L)[:L]
                d = d - np.bincount(we[m] + 1 - w0, minlength=L)[:L]
                cov[h] = (w0, np.cumsum(d)[:L - 1])
    return {"c0": c0, "cov": cov, "labs": labs}


def _rr_worker(r):
    P = _G
    parts = []
    for path, offs in P["mm_producers"]:
        a, b = int(offs[r]), int(offs[r + 1])
        if b > a:
            parts.append(np.load(path, mmap_mode="r")[a:b])
    lo = r * REGION_W
    hi = min(lo + REGION_W, P["mm_total"])
    flags = np.zeros(hi - lo, dtype=np.uint8)
    if not parts:
        return r, flags
    keys = np.sort(np.concatenate(parts))
    gw = keys >> MM_BITS
    vals = (keys & MM_MASK) - MM_OFF
    u, starts, counts = np.unique(gw, return_index=True, return_counts=True)
    ends = starts + counts
    bg = P["bg_mm"]
    ok = counts >= 5
    u, starts, counts, ends = u[ok], starts[ok], counts[ok], ends[ok]
    if len(u):
        fourth = vals[ends - 4]
        med_thr = np.maximum(5, np.minimum((counts * 0.1).astype(np.int64), 15))
        kth = vals[starts + med_thr - 1]
        mx = vals[ends - 1]
        f = (fourth >= bg) & (kth < bg) & (mx > bg)
        flags[(u[f] - lo).astype(np.int64)] = 1
    return r, flags


def _count_neighbourhoods(U, C, core, eb):
    """N(c): sum of C over the distinct +/-1-bin neighbour codes of each code in U[core]."""
    u = U[core]
    lo, hi = u >> eb, u & ((1 << eb) - 1)
    nu = len(U)
    nbs = []
    for da in (-2, 0, 2):
        for db in (-2, 0, 2):
            a, b = lo + da, hi + db
            nb = (np.minimum(a, b) << eb) | np.maximum(a, b)
            nb[(a < 0) | (b < 0)] = -1
            for prev in nbs:
                nb[nb == prev] = -1
            nbs.append(nb)
    N = np.zeros(len(u), dtype=np.int64)
    for nb in nbs:
        idx = np.searchsorted(U, nb)
        np.minimum(idx, nu - 1, out=idx)
        hit = (nb >= 0) & (U[idx] == nb)
        N += np.where(hit, C[idx], 0)
    return N


def _rc_worker(p):
    """Recurrence counts for codes whose first breakend lies in partition p."""
    P = _G
    eb = P["ebits"]
    b0, b1 = P["rc_bounds"][p], P["rc_bounds"][p + 1]
    klo = max(b0 - 2, 0) << eb
    khi = (b1 + 2) << eb
    parts = []
    for path in P["rc_files"]:
        a = np.load(path, mmap_mode="r")
        i, j = np.searchsorted(a, [klo, khi])
        if j > i:
            parts.append(np.array(a[i:j]))
    if not parts:
        return p, np.zeros(0, np.int64), np.zeros(0, np.int64)
    arr = np.concatenate(parts)
    del parts
    arr.sort()
    U, C = np.unique(arr, return_counts=True)
    del arr
    lo = U >> eb
    core = (lo >= b0) & (lo < b1)
    N = _count_neighbourhoods(U, C.astype(np.int64), core, eb)
    m = N >= 2
    return p, U[core][m], N[m]


def _recurrent_codes(code_files, ebits, threads):
    total = sum(n for _, n, _ in code_files)
    if total == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), 0
    sample = np.sort(np.concatenate([smp for _, _, smp in code_files]))
    nparts = int(max(threads * 4, -(-total // 4000000)))
    qs = np.unique(sample[np.linspace(0, len(sample) - 1, nparts + 1).astype(np.int64)[1:-1]])
    bounds = [0] + [int(x) for x in qs] + [1 << ebits]
    _G.clear()
    _G.update(ebits=ebits, rc_bounds=bounds, rc_files=[f for f, _, _ in code_files])
    outU, outN = {}, {}
    with _fork_pool(threads) as pool:
        for p, u, n in pool.imap_unordered(_rc_worker, range(len(bounds) - 1), chunksize=1):
            outU[p], outN[p] = u, n
    keys = sorted(outU)
    U = np.concatenate([outU[k] for k in keys]) if keys else np.zeros(0, np.int64)
    N = np.concatenate([outN[k] for k in keys]) if keys else np.zeros(0, np.int64)
    for f, _, _ in code_files:
        os.remove(f)
    return U, N, total


# ------------------------------------------------------------- patched Severus stages

def _setup(ref_lengths, args):
    lm = CTX["lm"]
    refnames = list(ref_lengths.keys())
    reflen = [int(ref_lengths[r]) for r in refnames]
    mm_nwin = [L // COV_WINDOW_MM + 2 for L in reflen]
    mm_off = np.concatenate([[0], np.cumsum(mm_nwin)]).astype(np.int64)
    n_cov = [L // COV_WINDOW + 1 for L in reflen]
    cov_off = np.concatenate([[0], np.cumsum(n_cov)]).astype(np.int64)
    parent = lm.tmp_dir or args.out_dir
    os.makedirs(parent, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="severus_lowmem_tmp_", dir=parent)
    W = lm.prescreen_window
    bin_bits = int(max(reflen) // W + 2).bit_length()
    tid_bits = int(len(refnames) + 1).bit_length()
    ebits = tid_bits + bin_bits + 1
    if 2 * ebits > 62:
        raise SystemExit(f"--prescreen-window {W} too small for this reference")
    CTX.update(
        refnames=refnames, reflen=reflen, mm_nwin=mm_nwin, mm_off=mm_off[:-1].tolist(),
        mm_total=int(mm_off[-1]), n_regions=int(-(-int(mm_off[-1]) // REGION_W)),
        cov_off=cov_off[:-1].tolist(), n_cov=int(cov_off[-1]), n_cov_list=n_cov,
        tmp=tmp, bin_bits=bin_bits, ebits=ebits, tidmap={r: i for i, r in enumerate(refnames)},
        bams=[], genome_ids=[], mm_producers=[], codes=[], bucket_files=[],
        cov=defaultdict(lambda: np.zeros(int(cov_off[-1]), dtype=np.int64)),
        gs_six=len(inspect.signature(sv_bp.get_segment).parameters) >= 6,
        lr_seven=len(inspect.signature(sv_bp.label_reads).parameters) >= 7,
        vntr_list=sv_vntr.read_vntr_file(args.vntr_file) if args.vntr_file else None,
        t0=time.time())
    prescreen = not lm.no_prescreen and lm.prescreen_min_reads > 1
    CTX["prescreen"] = prescreen
    logger.info(f"severus_lowmem {__version__} (Severus {SEVERUS_VERSION}); scratch: {tmp}")
    if prescreen:
        logger.info(f"\tPrescreen on: split reads need a junction shared by >= "
                    f"{lm.prescreen_min_reads} reads within +/-1 x {W} bp bins")
        if args.single_bp:
            logger.warning("\t--single-bp with prescreen: hanging breakpoints whose partner "
                           "ends are dispersed are withheld; use --no-prescreen to keep them")
        if args.write_alignments:
            logger.warning("\t--write-alignments with prescreen: withheld split reads are "
                           "not written")
        if getattr(args, "whitelist", None):
            logger.warning("\t--whitelist with prescreen: junctions seen in a single read are "
                           "withheld, including inside whitelisted regions")
        if args.bp_min_support < lm.prescreen_min_reads:
            logger.warning(f"\t--min-support {args.bp_min_support} is below "
                           f"--prescreen-min-reads {lm.prescreen_min_reads}")
    else:
        logger.info("\tPrescreen off: all reads are passed to the breakpoint caller")
    if SEVERUS_VERSION not in TESTED_VERSIONS:
        logger.warning(f"\tseverus_lowmem was tested with Severus {', '.join(TESTED_VERSIONS)}")


def _base_params(args):
    return dict(sv_size=args.sv_size, use_supp=args.use_supplementary_tag,
                ignore_hp=getattr(args, "ignore_hp", False), tmp=CTX["tmp"],
                refnames=CTX["refnames"], reflen=CTX["reflen"], mm_off=CTX["mm_off"],
                mm_nwin=CTX["mm_nwin"], n_regions=CTX["n_regions"], tidmap=CTX["tidmap"],
                cov_off=CTX["cov_off"], n_cov=CTX["n_cov"], gs_six=CTX["gs_six"],
                lr_seven=CTX["lr_seven"], prescreen=CTX["prescreen"],
                pw=CTX["lm"].prescreen_window, pmin=CTX["lm"].prescreen_min_reads,
                bin_bits=CTX["bin_bits"], ebits=CTX["ebits"], vntr_list=CTX["vntr_list"],
                min_mapq=args.min_mapping_quality, multisample=args.multisample,
                min_len=getattr(args, "min_len", -1))


def _log_progress(msg):
    logger.info(f"\t[lowmem] {msg} ({time.time() - CTX['t0']:.0f} s, "
                f"driver peak RSS {_rss_gb():.1f} GB)")


def lowmem_get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id,
                                  coverage_histograms, mismatch_histograms, n90ls, bg_mmls,
                                  read_qual, read_qual_len, args):
    if "refnames" not in CTX:
        _setup(ref_lengths, args)
    lm = CTX["lm"]
    bam_idx = len(CTX["bams"])
    CTX["bams"].append(bam_file)
    CTX["genome_ids"].append(genome_id)

    with pysam.AlignmentFile(bam_file, "rb") as fh:
        all_reference_ids = list(fh.references)
        try:
            mapped = fh.mapped
        except ValueError:
            mapped = 0
    if not mapped:
        mapped = os.path.getsize(bam_file) // 300   # conservative records estimate
    nb = lm.buckets or int(min(8192, max(4, -(-mapped // TARGET_BUCKET_ROWS))))
    fetch_list = []
    for j, ctg in enumerate(all_reference_ids):
        ctg_len = ref_lengths[ctg]
        for i in range(0, max(ctg_len // CHUNK_SIZE, 1)):
            reg_start = i * CHUNK_SIZE
            reg_end = (i + 1) * CHUNK_SIZE
            if ctg_len - reg_end < CHUNK_SIZE:
                reg_end = ctg_len
            fetch_list.append((j, ctg, reg_start, reg_end))
    cg0 = CTX.get("next_cg", 0)
    CTX["next_cg"] = cg0 + len(fetch_list)
    tasks = [(bam_idx, cg0 + i, j, ctg, s, e) for i, (j, ctg, s, e) in enumerate(fetch_list)]

    _G.clear()
    _G.update(_base_params(args))
    _G.update(bam_idx=bam_idx, bam_path=bam_file, genome_id=genome_id, nbuckets=nb)
    res = []
    with _fork_pool(args.threads) as pool:
        for r in pool.imap_unordered(_p1_worker, tasks, chunksize=1):
            res.append(r)
    res.sort(key=lambda r: r["cg"])
    n_rows = sum(r["n_rows"] for r in res)
    n_seg = sum(r["n_seg"] for r in res)
    _log_progress(f"parsed {len(res)} chunks: {n_rows:,} compact rows, {n_seg:,} segments, "
                  f"{nb} buckets; worker peak RSS {max(r['rss'] for r in res):.2f} GB")

    bgmm, st_rl, st_al, st_err, st_mm = {}, {}, {}, {}, {}
    for r in res:
        _merge_hist(bgmm, *r["bgmm"])
        _merge_hist(st_rl, *r["st_rl"])
        _merge_hist(st_al, *r["st_al"])
        _merge_hist(st_err, *r["st_err"])
        _merge_hist(st_mm, *r["st_mm"])
        if r["mm_path"]:
            CTX["mm_producers"].append((r["mm_path"], r["mm_offs"]))

    # bucket pass 1
    _G["chunk_files"] = [{k: r[k] for k in ("seg_path", "seg_off", "obj_off", "nam_off")}
                         for r in res]
    b1 = []
    with _fork_pool(args.threads) as pool:
        for r in pool.imap_unordered(_b1_worker, range(nb), chunksize=1):
            b1.append(r)
    b1.sort(key=lambda r: r["k"])
    n_reads = 0
    for r in b1:
        n_reads += r["n_reads"]
        _merge_hist(st_rl, *r["rl"])
        _merge_hist(st_err, *r["err"])
        _merge_hist(st_mm, *r["mm"])
        _merge_hist(st_al, *r["al"])
        CTX["mm_producers"].append((r["mm_path"], r["mm_offs"]))
        if r["codes_path"]:
            CTX["codes"].append((r["codes_path"], r["n_codes"], r["codes_sample"]))
    CTX["bucket_files"].append((bam_idx, nb, _G["chunk_files"]))
    _log_progress(f"regrouped {n_reads:,} reads with segments; worker peak RSS "
                  f"{max(r['rss'] for r in b1):.2f} GB")

    # background mismatch rate (background_mm_rat)
    QT = 0.95 if not args.multisample else 0.975
    bg_mm = _quantile_from_hist(*_hist_arrays(bgmm), QT)

    # read statistics (get_read_statistics)
    v_rl, c_rl = _hist_arrays(st_rl)
    v_al, c_al = _hist_arrays(st_al)
    reads_total_len = np.int64((v_rl * c_rl).sum())
    alignment_total_len = np.int64((v_al * c_al).sum())
    aln_rate = alignment_total_len / reads_total_len
    reads_n50 = _calc_nx_from_hist(v_rl, c_rl, reads_total_len, 0.50)
    reads_n90 = _calc_nx_from_hist(v_rl, c_rl, reads_total_len, 0.90)
    aln_n50 = _calc_nx_from_hist(v_al, c_al, alignment_total_len, 0.50)
    aln_n90 = _calc_nx_from_hist(v_al, c_al, alignment_total_len, 0.90)
    he, ce_ = _hist_arrays(st_err)
    hm, cm = _hist_arrays(st_mm)
    K_MM = sv_bp.K_MM
    e25, e50, e75 = (_quantile_from_hist(he, ce_, q) / K_MM for q in (0.25, 0.50, 0.75))
    m25, m50, m75 = (_quantile_from_hist(hm, cm, q) / K_MM for q in (0.25, 0.50, 0.75))
    logger.info(f"\tTotal read length: {reads_total_len}")
    logger.info(f"\tTotal aligned length: {alignment_total_len} ({aln_rate:.2f})")
    logger.info(f"\tRead N50 / N90: {reads_n50} / {reads_n90}")
    logger.info(f"\tAlignments N50 / N90: {aln_n50} / {aln_n90}")
    logger.info(f"\tRead error rate (Q25 / Q50 / Q75): {e25:.4f} / {e50:.4f} / {e75:.4f}")
    logger.info(f"\tRead mismatch rate (Q25 / Q50 / Q75): {m25:.4f} / {m50:.4f} / {m75:.4f}")
    n90 = min(reads_n90, args.min_aligned_length) if not args.multisample else args.min_aligned_length

    # coverage and read_qual from compact rows (update_cov_hist)
    _G.update(bg_mm_bam=bg_mm, n90_bam=n90)
    tasks = [(r["rows_path"], r["rows_chr"], r["cg"]) for r in res if r["rows_path"]]
    labs = {}
    cov = CTX["cov"]
    with _fork_pool(args.threads) as pool:
        for r in pool.imap_unordered(_ric_worker, tasks, chunksize=1):
            base = CTX["cov_off"][r["c0"]]
            for h, (w0, arr) in r["cov"].items():
                cov[(bam_idx, h)][base + w0: base + w0 + len(arr)] += arr
            for lab, (cnt, ln, first) in r["labs"].items():
                if lab in labs:
                    labs[lab] = (labs[lab][0] + cnt, labs[lab][1] + ln, min(labs[lab][2], first))
                else:
                    labs[lab] = (cnt, ln, first)
    for lab in sorted(labs, key=lambda x: labs[x][2]):
        read_qual[lab] += labs[lab][0]
        read_qual_len[lab] += labs[lab][1]
    for t in tasks:
        os.remove(t[0])

    n90ls.append(n90)
    bg_mmls.append(bg_mm)
    _log_progress(f"scratch in use {_dir_gb(CTX['tmp']):.1f} GB")
    return []


def lowmem_update_segments_by_read(segments_by_read, mismatch_histograms, bg_mm, ref_lengths,
                                   read_qual, read_qual_len, args):
    lm = CTX["lm"]
    bg = float(np.median(bg_mm))
    logger.info("Annotating reads")

    # junction recurrence (prescreen)
    rc_codes = rc_n = np.zeros(0, np.int64)
    if CTX["prescreen"]:
        rc_codes, rc_n, n_obs = _recurrent_codes(CTX.pop("codes"), CTX["ebits"], args.threads)
        _log_progress(f"prescreen: {n_obs:,} read-junctions, {len(rc_codes):,} junction bins "
                      f"shared by >= 2")

    # mismatch-rate windows (background_mm_hist), pooled over all BAMs
    _G.clear()
    _G.update(mm_producers=CTX["mm_producers"], mm_total=CTX["mm_total"], bg_mm=bg)
    flags = np.zeros(CTX["mm_total"], dtype=np.uint8)
    with _fork_pool(args.threads) as pool:
        for r, f in pool.imap_unordered(_rr_worker, range(CTX["n_regions"]), chunksize=1):
            flags[r * REGION_W: r * REGION_W + len(f)] = f
    for path, _ in CTX["mm_producers"]:
        os.remove(path)
    mm_hist_high = {}
    for c, name in enumerate(CTX["refnames"]):
        o = CTX["mm_off"][c]
        mm_hist_high[name] = flags[o:o + CTX["mm_nwin"][c]].tolist()
    del flags
    if args.write_segdups_out:
        sv_bp.extract_segdups(mm_hist_high, args.write_segdups_out)
    _log_progress(f"mismatch windows done (bg mismatch {bg})")

    # bucket pass 2 for each BAM
    seg_labs = {}
    parts = []
    tot = defaultdict(int)
    for bam_idx, nb, chunk_files in CTX["bucket_files"]:
        _G.clear()
        _G.update(_base_params(args))
        _G.update(bam_idx=bam_idx, genome_id=CTX["genome_ids"][bam_idx], nbuckets=nb,
                  chunk_files=chunk_files, bg_mm=bg, mm_hist_high=mm_hist_high,
                  min_aligned_length=args.min_aligned_length, rc_codes=rc_codes, rc_n=rc_n)
        segdiff = {}
        with _fork_pool(args.threads) as pool:
            for r in pool.imap_unordered(_b2_worker, range(nb), chunksize=1):
                for h, d in r["cov"].items():
                    if h not in segdiff:
                        segdiff[h] = np.zeros(CTX["n_cov"] + 1, dtype=np.int64)
                    segdiff[h] += d
                for lab, cnt in r["rq_n"].items():
                    first = r["rq_first"][lab]
                    if lab in seg_labs:
                        a = seg_labs[lab]
                        seg_labs[lab] = (a[0] + cnt, a[1] + r["rq_len"][lab], min(a[2], first))
                    else:
                        seg_labs[lab] = (cnt, r["rq_len"][lab], first)
                parts.append(r["path"])
                for x in ("n_full", "n_res", "n_withheld"):
                    tot[x] += r[x]
                tot["rss"] = max(tot["rss"], r["rss"])
        for ch in chunk_files:
            os.remove(ch["seg_path"])
        for h, d in segdiff.items():
            CTX["cov"][(bam_idx, h)] += np.cumsum(d)[:CTX["n_cov"]]
        del segdiff

    # read_qual.txt (write_readqual)
    for lab in sorted(seg_labs, key=lambda x: seg_labs[x][2]):
        read_qual[lab] += seg_labs[lab][0]
        read_qual_len[lab] += seg_labs[lab][1]
    with open(args.outpath_readqual, "w") as f:
        f.write('Number of segments:')
        f.writelines('{}\t{}\n'.format(k, v) for k, v in read_qual.items())
        f.write('Total length of segments:')
        f.writelines('\n{}\t{}'.format(k, v) for k, v in read_qual_len.items())

    # reads handed to the breakpoint caller, in Severus order
    items = []
    for p in parts:
        with open(p, "rb") as fh:
            items.extend(pickle.load(fh))
        os.remove(p)
    items.sort(key=lambda x: x[0])
    segments_by_read[:] = [r for _, r in items]
    del items
    if CTX["prescreen"]:
        _log_progress(f"{tot['n_full']:,} reads passed to the caller; junctions withheld for "
                      f"{tot['n_withheld']:,} split reads ({tot['n_res']:,} of them keep "
                      f"insertion/clipped segments); worker peak RSS {tot['rss']:.2f} GB")
    else:
        _log_progress(f"{tot['n_full']:,} reads passed to the caller; worker peak RSS "
                      f"{tot['rss']:.2f} GB")


_ORIG_update_coverage_hist = sv_bp.update_coverage_hist
_ORIG_extract_clipped_end = sv_bp.extract_clipped_end


def lowmem_update_coverage_hist(coverage_histograms, genome_ids, ref_lengths, segments_by_read,
                                control_genomes, target_genomes, loh_out):
    gid = CTX["genome_ids"]
    for (bam_idx, h), arr in CTX["cov"].items():
        for c, name in enumerate(CTX["refnames"]):
            key = (gid[bam_idx], h, name)
            o = CTX["cov_off"][c]
            add = arr[o:o + CTX["n_cov_list"][c]]
            if add.any():
                lst = coverage_histograms[key]
                coverage_histograms[key] = (np.asarray(lst, dtype=np.int64) + add).tolist()
    CTX["cov"].clear()
    return _ORIG_update_coverage_hist(coverage_histograms, genome_ids, ref_lengths, [],
                                      control_genomes, target_genomes, loh_out)


def lowmem_extract_clipped_end(segments_by_read):
    return _ORIG_extract_clipped_end([r for r in segments_by_read
                                      if not isinstance(r, ResidualRead)])


# --------------------------------------------------------------------------- main

def _lowmem_parser():
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--tmp-dir", default=None)
    p.add_argument("--no-prescreen", action="store_true")
    p.add_argument("--prescreen-window", type=int, default=500)
    p.add_argument("--prescreen-min-reads", type=int, default=2)
    p.add_argument("--buckets", type=int, default=0)
    p.add_argument("--keep-tmp", action="store_true")
    p.add_argument("--lowmem-help", action="store_true")
    return p


def main():
    lm, rest = _lowmem_parser().parse_known_args(sys.argv[1:])
    if lm.lowmem_help:
        print(__doc__)
        return 0
    if lm.prescreen_window < 1:
        sys.exit("--prescreen-window must be >= 1")
    if EXTRA_SLOTS:
        print(f"WARNING: unrecognised ReadSegment fields {list(EXTRA_SLOTS)}; they are "
              "carried as Python objects", file=sys.stderr)
    CTX.clear()
    CTX["lm"] = lm
    sv_main.get_all_reads_parallel = lowmem_get_all_reads_parallel
    sv_main.update_segments_by_read = lowmem_update_segments_by_read
    sv_main.update_coverage_hist = lowmem_update_coverage_hist
    sv_bf.extract_clipped_end = lowmem_extract_clipped_end
    sys.argv = [sys.argv[0]] + rest
    try:
        rc = sv_main.main()
    finally:
        tmp = CTX.get("tmp")
        if tmp and os.path.isdir(tmp) and not lm.keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    if CTX.get("t0"):
        logger.info(f"severus_lowmem finished in {time.time() - CTX['t0']:.0f} s; driver peak "
                    f"RSS {_rss_gb():.1f} GB")
    return rc


if __name__ == "__main__":
    sys.exit(main())
