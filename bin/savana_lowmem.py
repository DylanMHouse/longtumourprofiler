#!/usr/bin/env python3
"""
savana_lowmem.py - memory-bounded driver for SAVANA 1.3.7 `run`

Why this exists
---------------
Stock `savana run` holds the WHOLE GENOME's PotentialBreakpoint objects in the
parent process before clustering starts, then forks workers that copy-on-write
inherit that set for clustering/calling/annotation. On high-coverage ONT tumour
BAMs the accumulated set can exceed several hundred GB, so the job OOMs no
matter how few threads you use.

This driver keeps SAVANA's algorithms EXACTLY as they are (it imports and calls
the same functions) but restructures the orchestration into two passes:

  Pass A (scan):   chunked get_potential_breakpoints over the genome. Each
                   worker pickles its breakpoints STRAIGHT TO DISK, bucketed by
                   start-contig, and sends back only file paths. The parent
                   never accumulates breakpoints.

  Pass B (solve):  for each contig in turn - load only that contig's spilled
                   buckets, cluster -> call -> annotate depth, write per-contig
                   VCF/BEDPE/TSV/FASTA fragments, then free everything before
                   moving to the next contig.

Peak parent memory therefore becomes ~ (one contig's breakpoints) + (coverage
arrays), instead of (whole genome's breakpoints). chr1 is the worst case, so
expect roughly a 10-12x reduction versus stock on a human genome.

Translocations are preserved. The full contig list is always passed as
`contig_order`, so breakpoints whose mate is on another chromosome are NOT
dropped (that is what would happen if you naively split into per-chromosome
`savana run` jobs with --contigs chrN). Breakpoints generated while scanning
chrA but keyed to chrB are spilled into chrB's bucket and picked up when chrB
is solved in Pass B.

Usage
-----
Drop into your conda env, e.g.:
    $CONDA_PREFIX/lib/python3.11/site-packages/savana/savana_lowmem.py
or anywhere on PATH. Then run it INSTEAD of `savana run`:

    python savana_lowmem.py \
        --tumour DET011.somatic.forSV.bam \
        --normal DET011.germline.forSV.bam \
        --ref GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        --outdir savana_out --sample DET011 \
        --contigs contigs.chr.hg38.txt \
        --length 100 --threads 8 --coverage_binsize 100 \
        --chunksize 10000000 --spill_dir $TMPDIR/savana_spill --debug

Output is identical in form to `savana run`:
    <sample>.sv_breakpoints.vcf   (bcftools-sorted)
    <sample>.sv_breakpoints.bedpe
    <sample>.sv_breakpoints_read_support.tsv
    <sample>.inserted_sequences.fa
Feed the VCF to `savana classify` exactly as before.

NOTE: --spill_dir should point at fast node-local scratch (e.g. $TMPDIR on
BEAR). It holds the pickled intermediate breakpoints; budget roughly the same
order as the peak RAM stock SAVANA would have needed. It is cleaned up on
success unless --keep_spill is passed.
"""

import argparse
import gc
import glob
import os
import pickle
import re
import sys
import time
from math import ceil, floor
from multiprocessing import Pipe, Process, sharedctypes

# ---------------------------------------------------------------------------
# IMPORTANT: this file is often dropped INSIDE the installed savana package
# directory. When python runs a script, it puts that script's directory first
# on sys.path - which makes `import savana.helper` resolve to the sibling
# module savana.py instead of the savana PACKAGE, giving:
#     ModuleNotFoundError: No module named 'savana.run'; 'savana' is not a package
# Removing our own directory from sys.path before importing fixes it wherever
# the script is placed.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_REAL_HERE = os.path.realpath(_HERE)
# realpath, not abspath: on systems where the env is reached through a symlink
# (e.g. /rds/homes/<user>/<project> -> /rds/projects/<project>) the SAME
# directory appears in sys.path under a different spelling, and a plain string
# comparison leaves it in place.
sys.path[:] = [p for p in sys.path
               if os.path.realpath(p or os.getcwd()) != _REAL_HERE]
# make sure the directory CONTAINING the package (site-packages) is importable
_PARENT = os.path.dirname(_REAL_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)
# if something already bound the name `savana` to the non-package module, drop it
_stale = sys.modules.get('savana')
if _stale is not None and not hasattr(_stale, '__path__'):
    del sys.modules['savana']

import pysam
import pysam.bcftools as bcftools

import savana.helper as helper
from savana.breakpoints import get_potential_breakpoints, call_breakpoints, compute_depth
from savana.clusters import cluster_breakpoints


# --------------------------------------------------------------------------
# Pass A: scan + spill
# --------------------------------------------------------------------------

def scan_worker(task, cov_arrays_for_contig, conn):
    """Run get_potential_breakpoints for one chunk, spill results to disk.

    Sends back {start_contig: path} instead of the objects themselves. This is
    the single most important difference from stock SAVANA: the multiprocessing
    pipe never carries the breakpoint graph, so there is no pickle-buffer spike
    in the worker and no accumulation in the parent.
    """
    written = {}
    try:
        potential_breakpoints = get_potential_breakpoints(
            task['aln_file'], task['is_cram'], task['ref'], task['length'],
            task['mapq'], task['label'], task['contig_order'], task['contig'],
            task['start_pos'], task['end_pos'], task['coverage_binsize'],
            cov_arrays_for_contig, task['keep_inv_artefact'],
            task['inv_artefact_distance'], task['single_bnd'],
            task['single_bnd_min_length'], task['single_bnd_max_mapq']
        )
        for start_contig, bps in potential_breakpoints.items():
            if not bps:
                continue
            # one file per (target contig, task) so concurrent workers never
            # write the same path
            path = os.path.join(
                task['spill_dir'],
                f"pb.{start_contig}.{task['task_id']}.pkl"
            )
            with open(path, 'wb') as fh:
                pickle.dump(bps, fh, protocol=pickle.HIGHEST_PROTOCOL)
            written[start_contig] = path
        del potential_breakpoints
        gc.collect()
    except Exception as exc:  # noqa: BLE001 - report and let parent decide
        conn.send({'__error__': f"{task['contig']}:{task['start_pos']}-{task['end_pos']} {exc}"})
        conn.close()
        return
    conn.send(written)
    conn.close()


def build_scan_tasks(aln_files, args, contigs, contig_lengths):
    """Chunk the genome exactly the way stock SAVANA does."""
    tasks = []
    task_id = 0
    for label, aln_file in aln_files.items():
        for contig, contig_length in contig_lengths.items():
            if contig not in contigs:
                continue
            starts = []
            if contig_length > args.chunksize:
                n_intervals = floor(contig_length / args.chunksize) + 1
                pos = 0
                for _ in range(n_intervals):
                    end = min(pos + args.chunksize, contig_length)
                    starts.append((pos, end))
                    pos = end + 1
            else:
                starts.append((0, contig_length))
            for start_pos, end_pos in starts:
                if start_pos >= contig_length:
                    continue
                tasks.append({
                    'aln_file': aln_file.filename,
                    'label': label,
                    'is_cram': args.is_cram,
                    'ref': args.ref,
                    'length': args.length,
                    'mapq': args.mapq,
                    'contig_order': contigs,   # FULL list - keeps translocations
                    'contig': contig,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'coverage_binsize': args.coverage_binsize,
                    'keep_inv_artefact': args.keep_inv_artefact,
                    'inv_artefact_distance': args.inv_artefact_distance,
                    'single_bnd': args.single_bnd,
                    'single_bnd_min_length': args.single_bnd_min_length,
                    'single_bnd_max_mapq': args.single_bnd_max_mapq,
                    'spill_dir': args.spill_dir,
                    'task_id': task_id,
                })
                task_id += 1
    return tasks


def run_scan(aln_files, args, contigs, contig_lengths, cov_arrays):
    """Drive Pass A with a bounded pool of forked workers."""
    tasks = build_scan_tasks(aln_files, args, contigs, contig_lengths)
    total = len(tasks)
    if args.debug:
        print(f'Pass A: {total} scan tasks over {args.threads} process(es)', flush=True)

    n_slots = max(1, min(args.threads, total))
    procs = [None] * n_slots
    pipes = [None] * n_slots
    done = 0
    errors = []

    while tasks or any(p is not None for p in procs):
        for i in range(n_slots):
            if procs[i] is None and tasks:
                task = tasks.pop(0)
                parent_conn, child_conn = Pipe(duplex=False)
                # pass only this contig's coverage arrays to the child
                cov_for_contig = cov_arrays[task['label']][task['contig']]
                proc = Process(target=scan_worker,
                               args=(task, cov_for_contig, child_conn))
                proc.start()
                child_conn.close()
                procs[i] = (proc, task)
                pipes[i] = parent_conn
        for i in range(n_slots):
            if procs[i] is None:
                continue
            proc, task = procs[i]
            conn = pipes[i]
            # NOTE: poll(0) not poll(0.05). A blocking poll costs its timeout
            # for every slot that is not ready, so a sweep of N slots could
            # stall dispatch by N*timeout seconds - which throttles throughput
            # badly once you raise --threads.
            if conn.poll(0):
                try:
                    result = conn.recv()
                except EOFError:
                    result = {'__error__': f"worker died on task {task['task_id']}"}
                proc.join()
                conn.close()
                procs[i] = None
                pipes[i] = None
                done += 1
                if '__error__' in result:
                    errors.append(result['__error__'])
                    print(f"ERROR in scan task: {result['__error__']}",
                          file=sys.stderr, flush=True)
                if args.debug and done % 25 == 0:
                    print(f'  ...{done}/{total} scan tasks complete', flush=True)
            elif not proc.is_alive():
                # died without sending (OOM-killed child, segfault, etc.)
                proc.join()
                conn.close()
                procs[i] = None
                pipes[i] = None
                done += 1
                msg = (f"task {task['task_id']} ({task['label']} "
                       f"{task['contig']}:{task['start_pos']}-{task['end_pos']}) "
                       f"exited without result (exitcode={proc.exitcode})")
                errors.append(msg)
                print(f'ERROR: {msg}', file=sys.stderr, flush=True)
        time.sleep(0.01)

    if errors:
        print(f'WARNING: {len(errors)} scan task(s) failed - output is INCOMPLETE',
              file=sys.stderr, flush=True)
    return errors


# --------------------------------------------------------------------------
# Pass B: solve one contig at a time
# --------------------------------------------------------------------------

def load_contig_breakpoints(spill_dir, contig):
    """Load and concatenate every spilled bucket for one contig."""
    breakpoints = []
    for path in sorted(glob.glob(os.path.join(spill_dir, f'pb.{contig}.*.pkl'))):
        with open(path, 'rb') as fh:
            breakpoints.extend(pickle.load(fh))
    return breakpoints


def solve_contig(contig, args, cov_arrays, ref_fasta, counter_start):
    """cluster -> call -> annotate for a single contig. Returns output strings."""
    breakpoints = load_contig_breakpoints(args.spill_dir, contig)
    if not breakpoints:
        return '', '', '', '', 0, None

    _, clusters = cluster_breakpoints(
        contig, breakpoints, args.min_reads_per_cluster,
        args.buffer, args.insertion_buffer
    )
    del breakpoints
    gc.collect()

    called, _ = call_breakpoints(
        contig, clusters, args.end_buffer, args.length,
        args.min_support, args.tumour_only
    )
    del clusters
    gc.collect()

    called = compute_depth(called, cov_arrays, args.coverage_binsize)

    bedpe_s, vcf_s, tsv_s, fa_s = '', '', '', ''
    count = counter_start
    example = called[0] if called else None
    for bp in called:
        bedpe_s += bp.as_bedpe(count)
        vcf_s += bp.as_vcf(ref_fasta)
        tsv_s += bp.as_read_support(count)
        if bp.breakpoint_notation in ['<INS>', '+', '-']:
            fa_s += bp.as_insertion_fasta(count)
        count += 1
    n = len(called)
    del called
    gc.collect()
    return bedpe_s, vcf_s, tsv_s, fa_s, n, example


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def solve_worker(contig, args, cov_arrays, conn):
    """Solve one contig in a child process, writing fragments to spill_dir.

    Fragments are written with contig-LOCAL ids starting at 0; the parent
    renumbers them in contig order when concatenating, so the final files are
    byte-identical to the serial path.
    """
    try:
        ref_fasta = pysam.FastaFile(args.ref)
        bedpe_s, vcf_s, tsv_s, fa_s, n, example = solve_contig(
            contig, args, cov_arrays, ref_fasta, 0)
        base = os.path.join(args.spill_dir, f'frag.{contig}')
        for suffix, payload in (('vcf', vcf_s), ('bedpe', bedpe_s),
                                ('tsv', tsv_s), ('fa', fa_s)):
            with open(f'{base}.{suffix}', 'w') as fh:
                fh.write(payload)
        # an example breakpoint is needed for the VCF header; send a tiny
        # pickle of just one
        example_path = None
        if example is not None:
            example_path = f'{base}.example.pkl'
            with open(example_path, 'wb') as fh:
                pickle.dump(example, fh, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:  # noqa: BLE001
        conn.send({'__error__': f'{contig}: {exc}'})
        conn.close()
        return
    conn.send({'contig': contig, 'n': n, 'base': base, 'example': example_path})
    conn.close()


_ID_RE = re.compile(r'ID_(\d+)')


def renumber(text, offset):
    """Shift every ID_<n> token by offset (ids appear in all four formats)."""
    if not offset:
        return text
    return _ID_RE.sub(lambda m: f'ID_{int(m.group(1)) + offset}', text)


def run_solve_parallel(contigs, args, cov_arrays):
    """Pass B across several contigs at once. Returns per-contig fragment info."""
    pending = list(contigs)
    n_slots = max(1, min(args.solve_threads, len(pending)))
    procs = [None] * n_slots
    pipes = [None] * n_slots
    info = {}
    errors = []
    while pending or any(p is not None for p in procs):
        for i in range(n_slots):
            if procs[i] is None and pending:
                contig = pending.pop(0)
                parent_conn, child_conn = Pipe(duplex=False)
                proc = Process(target=solve_worker,
                               args=(contig, args, cov_arrays, child_conn))
                proc.start()
                child_conn.close()
                procs[i] = (proc, contig)
                pipes[i] = parent_conn
        for i in range(n_slots):
            if procs[i] is None:
                continue
            proc, contig = procs[i]
            conn = pipes[i]
            if conn.poll(0):
                try:
                    res = conn.recv()
                except EOFError:
                    res = {'__error__': f'{contig}: worker died'}
                proc.join()
                conn.close()
                procs[i] = None
                pipes[i] = None
                if '__error__' in res:
                    errors.append(res['__error__'])
                    print(f"ERROR solving {res['__error__']}", file=sys.stderr, flush=True)
                else:
                    info[contig] = res
                    if args.debug:
                        print(f"  {contig}: {res['n']} breakpoints", flush=True)
            elif not proc.is_alive():
                proc.join()
                conn.close()
                procs[i] = None
                pipes[i] = None
                msg = f'{contig}: exited without result (exitcode={proc.exitcode})'
                errors.append(msg)
                print(f'ERROR solving {msg}', file=sys.stderr, flush=True)
        time.sleep(0.01)
    return info, errors


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Memory-bounded per-contig driver for SAVANA run',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--tumour', required=True)
    p.add_argument('--normal', required=False, default=None)
    p.add_argument('--ref', required=True)
    p.add_argument('--ref_index', default=None,
                   help='defaults to <ref>.fai')
    p.add_argument('--outdir', required=True)
    p.add_argument('--sample', default=None,
                   help='defaults to tumour basename')
    p.add_argument('--contigs', default=None)
    p.add_argument('--spill_dir', default=None,
                   help='scratch for spilled breakpoints (default <outdir>/spill)')
    p.add_argument('--keep_spill', action='store_true')
    # algorithm params - defaults copied from savana 1.3.7 `run`
    p.add_argument('--length', type=int, default=30)
    p.add_argument('--mapq', type=int, default=5)
    p.add_argument('--buffer', type=int, default=10)
    p.add_argument('--insertion_buffer', type=int, default=250)
    p.add_argument('--end_buffer', type=int, default=50)
    p.add_argument('--min_support', type=int, default=3)
    p.add_argument('--min_reads_per_cluster', type=int, default=3)
    p.add_argument('--coverage_binsize', type=int, default=100)
    p.add_argument('--chunksize', type=int, default=10000000)
    p.add_argument('--threads', type=int, default=8,
                   help='processes for Pass A (BAM scanning)')
    p.add_argument('--solve_threads', type=int, default=1,
                   help='contigs solved concurrently in Pass B. Each holds one '
                        'contig in RAM, so peak = solve_threads x largest contig. '
                        '1 = lowest memory (default).')
    p.add_argument('--single_bnd', action='store_true')
    p.add_argument('--single_bnd_min_length', type=int, default=1000)
    p.add_argument('--single_bnd_max_mapq', type=int, default=5)
    p.add_argument('--keep_inv_artefact', action='store_true')
    p.add_argument('--inv_artefact_distance', type=int, default=200)
    p.add_argument('--debug', action='store_true')
    # Accepted for command-line compatibility with `savana run`, but NOT used
    # by SV breakpoint calling. In savana 1.3.7 --snp_vcf feeds allele_counter
    # and --blacklist feeds bin_generator/read_counter - both are copy-number
    # only. They have never had any effect on the breakpoint/SV path.
    p.add_argument('--snp_vcf', default=None, help=argparse.SUPPRESS)
    p.add_argument('--blacklist', default=None, help=argparse.SUPPRESS)
    p.add_argument('--g1000_vcf', default=None, help=argparse.SUPPRESS)
    args = p.parse_args(argv)

    for ignored in ('snp_vcf', 'blacklist', 'g1000_vcf'):
        if getattr(args, ignored):
            print(f'NOTE: --{ignored} is a copy-number option in savana and has '
                  f'no effect on SV calling; ignoring it.', file=sys.stderr)
            setattr(args, ignored, None)

    if args.ref_index is None:
        args.ref_index = args.ref + '.fai'
    if args.sample is None:
        args.sample = os.path.basename(args.tumour).rsplit('.bam', 1)[0].rsplit('.cram', 1)[0]
    if args.spill_dir is None:
        args.spill_dir = os.path.join(args.outdir, 'spill')
    args.is_cram = args.tumour.endswith('.cram')
    args.tumour_only = args.normal is None
    # stock savana_run does this clamp; replicate it so results match
    if args.min_reads_per_cluster > args.min_support:
        print(f'Lowering `min_reads_per_cluster` to {args.min_support} '
              f'to match `min_support`', file=sys.stderr)
        args.min_reads_per_cluster = args.min_support
    return args


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.spill_dir, exist_ok=True)

    if not os.path.exists(args.ref_index):
        sys.exit(f'ERROR: reference index not found: {args.ref_index}')

    aln_files = {'tumour': pysam.AlignmentFile(args.tumour, 'rc' if args.is_cram else 'rb')}
    if args.normal:
        aln_files['normal'] = pysam.AlignmentFile(
            args.normal, 'rc' if args.is_cram else 'rb')

    contigs = helper.get_contigs(args.contigs, args.ref_index)
    contig_lengths = helper.get_contig_lengths(args.ref_index)
    contigs = [c for c in contigs if c in contig_lengths]

    # coverage arrays: shared, genome-wide, cheap at binsize>=100
    approx_gb = sum(
        ceil(contig_lengths[c] / args.coverage_binsize) * 4 * 3
        for c in contigs) * len(aln_files) / 1e9
    print(f'Allocating coverage arrays (~{approx_gb:.2f} GB at '
          f'binsize {args.coverage_binsize})', flush=True)
    cov_arrays = {}
    for label in aln_files:
        for contig in contigs:
            n_bins = ceil(contig_lengths[contig] / args.coverage_binsize)
            cov_arrays.setdefault(label, {})[contig] = {
                hp: sharedctypes.RawArray('i', n_bins) for hp in [1, 2, None]
            }

    t0 = time.time()
    errors = run_scan(aln_files, args, contigs, contig_lengths, cov_arrays)
    print(f'Pass A (scan + spill) complete in {time.time() - t0:.1f}s', flush=True)

    # Pass B
    ref_fasta = pysam.FastaFile(args.ref)
    vcf_file = os.path.join(args.outdir, f'{args.sample}.sv_breakpoints.vcf')
    bedpe_file = os.path.join(args.outdir, f'{args.sample}.sv_breakpoints.bedpe')
    tsv_file = os.path.join(args.outdir, f'{args.sample}.sv_breakpoints_read_support.tsv')
    fa_file = os.path.join(args.outdir, f'{args.sample}.inserted_sequences.fa')

    t1 = time.time()
    total_bp = 0
    example_bp = None
    body_path = vcf_file + '.body.tmp'
    solve_errors = []

    if args.solve_threads > 1:
        # --- parallel Pass B ---
        info, solve_errors = run_solve_parallel(contigs, args, cov_arrays)
        with open(body_path, 'w') as vcf_body, \
             open(bedpe_file, 'w') as bedpe_out, \
             open(tsv_file, 'w') as tsv_out, \
             open(fa_file, 'w') as fa_out:
            tsv_out.write('VARIANT_ID\tTUMOUR_SUPPORTING_READS\tNORMAL_SUPPORTING_READS\n')
            for contig in contigs:          # contig order preserves numbering
                res = info.get(contig)
                if not res or not res['n']:
                    continue
                if example_bp is None and res['example']:
                    with open(res['example'], 'rb') as fh:
                        example_bp = pickle.load(fh)
                for suffix, handle in (('vcf', vcf_body), ('bedpe', bedpe_out),
                                       ('tsv', tsv_out), ('fa', fa_out)):
                    frag = f"{res['base']}.{suffix}"
                    if os.path.exists(frag):
                        with open(frag) as fh:
                            handle.write(renumber(fh.read(), total_bp))
                total_bp += res['n']
        for res in info.values():
            for suffix in ('vcf', 'bedpe', 'tsv', 'fa', 'example.pkl'):
                frag = f"{res['base']}.{suffix}"
                if os.path.exists(frag):
                    os.remove(frag)
    else:
        # --- serial Pass B (lowest memory) ---
        with open(body_path, 'w') as vcf_body, \
             open(bedpe_file, 'w') as bedpe_out, \
             open(tsv_file, 'w') as tsv_out, \
             open(fa_file, 'w') as fa_out:
            tsv_out.write('VARIANT_ID\tTUMOUR_SUPPORTING_READS\tNORMAL_SUPPORTING_READS\n')
            for contig in contigs:
                bedpe_s, vcf_s, tsv_s, fa_s, n, example = solve_contig(
                    contig, args, cov_arrays, ref_fasta, total_bp)
                if example is not None and example_bp is None:
                    example_bp = example
                vcf_body.write(vcf_s)
                bedpe_out.write(bedpe_s)
                tsv_out.write(tsv_s)
                fa_out.write(fa_s)
                total_bp += n
                if args.debug:
                    print(f'  {contig}: {n} breakpoints', flush=True)
                del bedpe_s, vcf_s, tsv_s, fa_s
                gc.collect()

    errors = errors + solve_errors

    print(f'Pass B (cluster/call/annotate) complete in {time.time() - t1:.1f}s', flush=True)
    print(f'{total_bp} annotated breakpoints to output', flush=True)

    # stitch header + body, then sort
    with open(vcf_file, 'w') as out:
        out.write(helper.generate_vcf_header(args, example_bp))
        with open(body_path) as body:
            for line in body:
                out.write(line)
    os.remove(body_path)
    bcftools.sort('-o', vcf_file, vcf_file, catch_stdout=False)

    if not args.keep_spill and not errors:
        for path in glob.glob(os.path.join(args.spill_dir, 'pb.*.pkl')):
            os.remove(path)
        try:
            os.rmdir(args.spill_dir)
        except OSError:
            pass
    elif errors:
        print(f'Spill files retained in {args.spill_dir} (errors occurred)',
              file=sys.stderr, flush=True)

    print(f'Wrote {vcf_file}', flush=True)
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
