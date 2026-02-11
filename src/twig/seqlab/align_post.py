#!/usr/bin/env python

"""Post-trimming actions

1. optionally subselect samples to keep/exclude
2. convert to AA
3. trim AA columns by gaps, and rows by min overlap, convert back to NT
X. optionally refine alignment
4. mask FS/STOP and write to NT & AA
5. write log

"""

from __future__ import annotations
from typing import Dict, List, Tuple
import sys
import re
import subprocess
from pathlib import Path
from loguru import logger
import toytree
import numpy as np
from twig.utils import set_log_level, TwigError

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")
BIN_TRIMAL = str(BIN / "trimal")
MISSING = set(["-", "N", "n", "?", "."])


def parse_fasta_to_dict(path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA starts with sequence before any header.")
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta_from_dict(path: Path, seqs: Dict[str, str]) -> None:
    """Write a FASTA dict to file."""
    with open(path, "w") as hout:
        for name, seq in seqs.items():
            hout.write(f">{name}\n{seq}\n")


def build_trimal_id_maps(names: List[str]) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return {original->temp} and {temp->original} ID maps for trimal-safe headers."""
    orig_to_tmp: Dict[str, str] = {}
    tmp_to_orig: Dict[str, str] = {}
    for idx, name in enumerate(names, start=1):
        tid = f"TRIMAL_{idx:06d}"
        orig_to_tmp[name] = tid
        tmp_to_orig[tid] = name
    return orig_to_tmp, tmp_to_orig


def write_renamed_fasta(src: Path, dst: Path, name_map: Dict[str, str]) -> None:
    """Write FASTA with headers replaced using name_map {old_name: new_name}."""
    seqs = parse_fasta_to_dict(src)
    out: Dict[str, str] = {}
    for name, seq in seqs.items():
        if name not in name_map:
            raise TwigError(f"missing header in remap table: {name}")
        out[name_map[name]] = seq
    write_fasta_from_dict(dst, out)


def restore_trimal_headers(seqs: Dict[str, str], tmp_to_orig: Dict[str, str]) -> Dict[str, str]:
    """Return sequence dict with temporary trimal IDs replaced by original headers."""
    restored: Dict[str, str] = {}
    for tname, seq in seqs.items():
        if tname not in tmp_to_orig:
            raise TwigError(f"unexpected trimal ID not found in mapping table: {tname}")
        oname = tmp_to_orig[tname]
        if oname in restored:
            raise TwigError(f"duplicate restored header after trimal remap: {oname}")
        restored[oname] = seq
    return restored


def get_alignment_stats(seqs: Dict[str, str], stage: str) -> Tuple[int, int]:
    """Return (nseqs, seqlen) for aligned fasta-like dict, else raise TwigError."""
    if not seqs:
        raise TwigError(f"no sequences remain after {stage}")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise TwigError(f"sequences are not aligned after {stage} (variable lengths)")
    return len(seqs), lengths.pop()


def run_checked(cmd: List[str], step: str):
    """Run subprocess command and raise TwigError with stderr/stdout on failure."""
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        msg = (proc.stderr or proc.stdout or "unknown error").strip()
        raise TwigError(f"{step} failed: {msg}")
    return proc


# def call_macse_refine_alignment(data: Path, outprefix: str, max_iter: int, verbose: bool):
#     """Run Alignment step with default settings"""
#     cmd = [
#         BIN_MACSE, "-prog", "refineAlignment",
#         "-align", str(data),
#         "-out_NT", f"{outprefix}.tmp.rmsa.nt.fa",
#         "-out_AA", f"{outprefix}.tmp.rmsa.aa.fa",
#         "-max_refine_iter", str(max_iter),
#     ]
#     logger.info("refining alignment")
#     logger.debug(f"[{data.name}] " + " ".join(cmd))
#     if verbose:
#         proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
#     else:
#         proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#     if proc.returncode:
#         raise Exception(proc.stderr)
#     return outprefix.with_suffix(outprefix.suffix + ".tmp.rmsa.nt.fa")


def filter_by_selection(fasta: Path, outprefix: Path, exclude: List[str], subsample: List[str], subsample_tree: Path) -> Tuple[Path, List[str]]:
    """Optionally filter aligned sequences by regex include/exclude or tree tips."""
    seqs = parse_fasta_to_dict(fasta)

    # raise an exception if not aligned
    lengths = [len(i) for i in seqs.values()]
    if not len(set(lengths)) == 1:
        raise TwigError("input contains sequences of variable lengths (not aligned)")

    # expand exclude list of names/globs to select all name matches
    matched = []
    names = list(seqs)
    if exclude:
        patterns = [re.compile(p) for p in exclude]
        matched = [s for s in names if any(r.search(s) for r in patterns)]
    if subsample:
        patterns = [re.compile(p) for p in subsample]
        matched = [s for s in names if any(r.search(s) for r in patterns)]
        matched = list(set(names) - set(matched))
    if subsample_tree:
        tree_tips = set(toytree.tree(subsample_tree).get_tip_labels())
        tips_missing_in_fasta = sorted(tree_tips - set(names))
        for tip in tips_missing_in_fasta:
            logger.warning(f"{tip} present in --subsample-tree but missing in fasta")
        matched = list(set(names) - tree_tips)
        for name in sorted(matched):
            logger.info(f"{name} excluded because absent from --subsample-tree tips")

    # keep only sequences that match the selection
    filtered = []
    keep = {}
    for name, seq in seqs.items():
        if name in matched:
            logger.debug(f"{name} excluded by user include/exclude args")
            filtered.append(name)
        else:
            keep[name] = seq

    # write output
    out = outprefix.with_suffix(outprefix.suffix + ".post")
    with open(out, 'w') as hout:
        for uname in keep:
            hout.write(f">{uname}\n{keep[uname]}\n")
    return out, filtered


def call_macse_trim_alignment(msa_nt: Path, edge_trim_window_size: int, edge_trim_min_coverage: float):
    """..."""
    out_nt = msa_nt.with_suffix(msa_nt.suffix + ".edge_trimmed")
    out_info = msa_nt.with_suffix(msa_nt.suffix + ".edge_trimmed.info")
    cmd = [
        BIN_MACSE, "-prog", "trimAlignment",
        "-align", str(msa_nt),
        "-respect_first_RF_ON",
        "-half_window_size", str(edge_trim_window_size),
        "-min_percent_NT_at_ends", str(edge_trim_min_coverage),
        "-out_trim_info", str(out_info),
        "-out_NT", str(out_nt),
    ]
    logger.debug(" ".join(cmd))
    run_checked(cmd, "MACSE trimAlignment")
    return out_nt


def call_macse_export_pre(msa_nt: Path):
    """..."""
    out_nt = msa_nt.with_suffix(msa_nt.suffix + ".null")
    out_aa = msa_nt.with_suffix(msa_nt.suffix + ".translated")
    cmd = [
        BIN_MACSE, "-prog", "exportAlignment",
        "-align", str(msa_nt),
        "-out_NT", str(out_nt),
        "-out_AA", str(out_aa),
    ]

    # write the NT and AA files
    run_checked(cmd, "MACSE exportAlignment (pre-trimal)")

    # remove gaps from NT file to emulate non-aligned
    seqs = parse_fasta_to_dict(out_nt)
    with open(out_nt, 'w') as hout:
        for uname in seqs:
            hout.write(f">{uname}\n{seqs[uname].replace('-', '')}\n")

    # return paths
    return out_aa, out_nt


def call_trimal(msa_aa: Path, iso_nt: Path, trimal_resoverlap: float, trimal_seqoverlap: float, trimal_algorithm: str):
    """..."""
    # trimal can truncate headers at ':' / ',' / whitespace; remap to temporary IDs first.
    aa_seqs = parse_fasta_to_dict(msa_aa)
    nt_seqs = parse_fasta_to_dict(iso_nt)
    aa_names = list(aa_seqs.keys())
    nt_names = list(nt_seqs.keys())
    aa_set = set(aa_names)
    nt_set = set(nt_names)
    if aa_set != nt_set:
        only_aa = sorted(aa_set - nt_set)
        only_nt = sorted(nt_set - aa_set)
        raise TwigError(
            "AA/NT headers differ before trimal remap: "
            f"AA-only={only_aa[:5]} NT-only={only_nt[:5]}"
        )

    orig_to_tmp, tmp_to_orig = build_trimal_id_maps(aa_names)
    safe_aa = msa_aa.with_suffix(msa_aa.suffix + ".safe")
    safe_nt = iso_nt.with_suffix(iso_nt.suffix + ".safe")
    write_renamed_fasta(msa_aa, safe_aa, orig_to_tmp)
    write_renamed_fasta(iso_nt, safe_nt, orig_to_tmp)
    logger.info(f"remapped {len(orig_to_tmp)} sequence headers to temporary trimal-safe IDs")

    out_nt = msa_aa.with_suffix(msa_aa.suffix + ".in_trimmed.backtranslated")
    cmd = [
        BIN_TRIMAL,
        "-in", str(safe_aa),
        "-out", str(out_nt),
        "-backtrans", str(safe_nt),
        "-fasta",
        "-ignorestopcodon",
        "-resoverlap", str(trimal_resoverlap),
        "-seqoverlap", str(int(100 * trimal_seqoverlap)),     # trimal wants this as an int for some reason
    ]
    if trimal_algorithm:
        cmd += [f"-{trimal_algorithm}"]
    logger.debug(" ".join(cmd))
    run_checked(cmd, "trimal backtranslation")

    # restore original headers and report filtered names using original IDs.
    post_tmp = parse_fasta_to_dict(out_nt)
    if not post_tmp:
        raise TwigError("no sequences remain after trimal trimming")
    pre_tmp_ids = set(orig_to_tmp.values())
    post_tmp_ids = set(post_tmp.keys())
    filtered_tmp = [tid for tid in orig_to_tmp.values() if tid not in post_tmp_ids]
    filtered = [tmp_to_orig[tid] for tid in filtered_tmp]

    restored = restore_trimal_headers(post_tmp, tmp_to_orig)
    write_fasta_from_dict(out_nt, restored)
    trim_len = len(list(restored.values())[0])
    return out_nt, filtered, trim_len


def call_macse_export_alignment_final(msa_nt: Path, outprefix: str):  #, codon_efs, codon_ifs, codon_fst, codon_ist):
    """..."""
    out_nt = outprefix.with_suffix(outprefix.suffix + ".nt.fa")
    out_aa = outprefix.with_suffix(outprefix.suffix + ".aa.fa")
    cmd = [
        BIN_MACSE, "-prog", "exportAlignment",
        "-align", str(msa_nt),
        "-codonForExternalFS",   "NNN", #str(codon_efs),
        "-codonForInternalFS",   "NNN", #str(codon_ifs),
        "-codonForFinalStop",    "NNN", #str(codon_fst),
        "-codonForInternalStop", "NNN", #str(codon_ist),
        "-out_NT", str(out_nt),
        "-out_AA", str(out_aa),
    ]
    logger.debug(" ".join(cmd))
    run_checked(cmd, "MACSE exportAlignment (final)")
    return out_nt, out_aa



###############################################################################
###############################################################################
###############################################################################


def present_matrix(seqs: Dict[str, str], missing=MISSING) -> Tuple[List[str], np.ndarray]:
    names = list(seqs.keys())
    L = len(seqs[names[0]])
    for n in names[1:]:
        if len(seqs[n]) != L:
            raise ValueError("All sequences must have the same alignment length.")
    arr = np.empty((len(names), L), dtype=np.bool_)
    for i, n in enumerate(names):
        s = seqs[n]
        arr[i] = np.fromiter((c not in missing for c in s), count=L, dtype=np.bool_)
    return names, arr


def overlap_matrix(present: np.ndarray) -> np.ndarray:
    x = present.astype(np.uint16, copy=False)
    return x @ x.T  # (n,n), diagonal = present counts


def prune_to_pairwise_min_overlap(names: List[str], overlap: np.ndarray, min_ov: int, present_counts: np.ndarray | None = None, min_keep: int = 4) -> Tuple[List[str], List[Tuple[str, int, str]], bool]:
    """
    Keep a set S such that overlap(i,j) >= min_ov for EVERY i!=j in S.
    Drop one sample at a time until satisfied, or fail (kept < min_keep).

    Returns:
      kept_names,
      removed_log: (removed_name, min_overlap_at_removal, partner_name_at_that_min),
      success: bool
    """
    n = len(names)
    if overlap.shape != (n, n):
        raise ValueError("overlap must be (n,n) matching names length.")

    active = np.ones(n, dtype=bool)
    removed: List[Tuple[str, int, str]] = []

    while active.sum() >= min_keep:
        idx = np.flatnonzero(active)
        sub = overlap[np.ix_(idx, idx)].astype(int, copy=False)

        k = sub.shape[0]
        if k < 2:
            break

        offdiag = ~np.eye(k, dtype=bool)

        # global minimum overlap among remaining pairs (excluding diagonal)
        min_val = sub[offdiag].min()
        if min_val >= min_ov:
            return [names[i] for i in idx], removed, True

        # one worst pair (i,j) achieving min_val
        wi, wj = np.argwhere((sub == min_val) & offdiag)[0]
        gi, gj = int(idx[wi]), int(idx[wj])

        # decide which endpoint to drop:
        # drop the one with smaller mean overlap to others (tie-break: fewer present sites)
        row_sums = sub.sum(axis=1) - np.diag(sub)         # exclude self
        means = row_sums / (k - 1)

        drop_local = wi if means[wi] < means[wj] else wj

        if means[wi] == means[wj] and present_counts is not None:
            pci, pcj = present_counts[gi], present_counts[gj]
            drop_local = wi if pci <= pcj else wj

        drop_global = int(idx[drop_local])
        partner_global = gj if drop_global == gi else gi

        removed.append((names[drop_global], int(min_val), names[partner_global]))
        active[drop_global] = False

    kept = [names[i] for i in np.flatnonzero(active)]
    return kept, removed, False


def filter_by_min_overlap(fasta: Path, min_ov: int, min_samples: int):
    seqs = parse_fasta_to_dict(fasta)
    if not seqs:
        raise TwigError("no sequences remain before min-overlap filtering")
    names, arr = present_matrix(seqs)
    ovarr = overlap_matrix(arr)
    kept, filtered, _ = prune_to_pairwise_min_overlap(names, ovarr, min_ov, arr.sum(axis=1), min_samples)

    # write to tmp file
    if filtered:
        out_nt = fasta.with_suffix(fasta.suffix + ".terrace.trimmed")
        with open(out_nt, 'w') as hout:
            for uname in seqs:
                if uname in kept:
                    hout.write(f">{uname}\n{seqs[uname]}\n")
        return out_nt, filtered, kept
    return fasta, filtered, kept


###############################################################################
###############################################################################
###############################################################################


def run_align_post(args):
    """..."""
    set_log_level(args.log_level)#, args.log_file)

    # check tools
    if not Path(BIN_MACSE).exists():
        raise TwigError(f"macse binary not found. Checked: {BIN_MACSE}")
    if (not args.skip_trimal) and (not Path(BIN_TRIMAL).exists()):
        raise TwigError(f"trimal binary not found. Checked: {BIN_TRIMAL}")

    # check infiles
    if not (args.input.exists() and args.input.is_file()):
        raise TwigError(f"{args.input} not found")
    if args.subsample_tree and not (args.subsample_tree.exists() and args.subsample_tree.is_file()):
        raise TwigError(f"{args.subsample_tree} not found")

    # only one or the other allowed
    nargs = len([i for i in [args.exclude, args.subsample, args.subsample_tree] if i])
    if nargs > 1:
        raise TwigError("choose one of --exclude, --subsample, or --subsample-tree")

    # ensure outpath and pdir exists
    if args.out_prefix.is_dir():
        raise TwigError("out-prefix must be a file path prefix, not a directory")
    args.out_prefix.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    result = args.out_prefix.with_suffix(args.out_prefix.suffix + ".nt.fa")
    if result.exists() and not args.force:
        logger.warning(f"[skipping] result {result} already exists. Use --force to overwrite")
        return
    ####################################################################

    try:
        # record starting size
        seqs0 = parse_fasta_to_dict(args.input)
        n0, s0 = get_alignment_stats(seqs0, "input parsing")
        logger.debug((n0, s0))

        # FILTER TO OPTIONALLY SUBSAMPLE -> tmp.msa.nt.fa
        msa_nt, f_user = filter_by_selection(args.input, args.out_prefix, args.exclude, args.subsample, args.subsample_tree)
        for sample in f_user:
            logger.info(f"{sample} filtered by user include/exclude args")
        n_sel, s_sel = get_alignment_stats(parse_fasta_to_dict(msa_nt), "selection filtering")

        # TRIM TO MIN-GAP EDGES -> tmp.msa.nt.fa
        msa_nt = call_macse_trim_alignment(msa_nt, args.edge_window_size, args.edge_min_coverage)

        # POST TRIMMING/FILTERING OF THE AA ALIGNMENT IN TRIMAL
        if args.skip_trimal:
            f_trim = []
            trim_len = s_sel
        else:
            # export as AA for trimming and run trimal to get NT result
            msa_aa, iso_nt = call_macse_export_pre(msa_nt)
            msa_nt, f_trim, trim_len = call_trimal(msa_aa, iso_nt, args.trimal_res_overlap, args.trimal_seq_overlap, args.trimal_algorithm)
            for sample in f_trim:
                logger.info(f"{sample} filtered by trimal seqoverlap")
            if trim_len < args.min_retained_nt_length:
                raise TwigError("trimmed length < min_retained_nt_length")

        seqsx = parse_fasta_to_dict(msa_nt)
        nx, sx = get_alignment_stats(seqsx, "AA/trimal filtering")
        if sx < args.min_retained_nt_length:
            raise TwigError("trimmed length < min_retained_nt_length")

        # APPLY MIN TERRACE OVERLAP
        msa_nt, f_terrace, kept = filter_by_min_overlap(msa_nt, args.min_pair_overlap, args.min_retained_seqs)
        if len(kept) < args.min_retained_seqs:
            raise TwigError(f"N passed samples ({len(kept)}) < min_retained_seqs")
        for sample in f_terrace:
            logger.info(f"{sample[0]} filtered by min overlap")
        n_terrace = len(kept)

        # [OPTIONAL] REFINE ALIGNMENT if any sequences were removed
        # if args.refine_alignment and  or (args.refine_alignment_if and filtered):
        #     data = call_macse_refine_alignment(data, args.out_prefix, args.max_iter_refine_alignment, args.verbose)

        # MASK FS,STOP and translate
        out_nt, out_aa = call_macse_export_alignment_final(msa_nt, args.out_prefix)  # args.codon_int_fs, args.codon_ext_fs, args.codon_final_stop, args.codon_int_stop)
        seqs1 = parse_fasta_to_dict(out_nt)
        n1, s1 = get_alignment_stats(seqs1, "final export")
        removed_user = max(0, n0 - n_sel)
        removed_trimal = max(0, n_sel - nx)
        removed_min_overlap = max(0, nx - n_terrace)
        removed_total = max(0, n0 - n1)
        if (not args.skip_trimal) and (len(f_trim) != removed_trimal):
            logger.warning(
                f"trimal removed-count mismatch (name-based={len(f_trim)} vs count-based={removed_trimal}); using count-based in summary"
            )
        logger.info(
            "sequence removals by stage: "
            f"user={removed_user}, trimal={removed_trimal}, min_overlap={removed_min_overlap}, total={removed_total}"
        )
        logger.info(f"sequence alignment trimmed from ({n0}, {s0} nt) to ({n1}, {s1} nt)")
        logger.info(f"post-trimmed NT alignment written to {out_nt}")
        logger.info(f"post-trimmed AA alignment written to {out_aa}")

    # clean up tmp files
    finally:
        suffices = [
            ".post",
            ".post.edge_trimmed",
            ".post.edge_trimmed.info",
            ".post.edge_trimmed.null",
            ".post.edge_trimmed.null.safe",
            ".post.edge_trimmed.translated",
            ".post.edge_trimmed.translated.safe",
            ".post.edge_trimmed.translated.in_trimmed.backtranslated",
            ".post.edge_trimmed.terrace.trimmed",
            ".post.edge_trimmed.translated.in_trimmed.backtranslated.terrace.trimmed",
        ]
        if not args.keep:
            for suffix in suffices:
                path = args.out_prefix.with_suffix(args.out_prefix.suffix + suffix)
                if path.exists():
                    path.unlink()


def main():
    from ..cli.subcommands import get_parser_align_post
    parser = get_parser_align_post()
    args = parser.parse_args()
    run_align_post(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
