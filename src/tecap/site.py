"""Read priming site and its label, on the definition of flight-polyafilter-bench.

terminal_a_run is a port of pafbench/bd67.py (sense_3p_clip, inward_pairs,
terminal_a_run) with the strand passed in: the run is searched at the
gene-strand 3' end of the BAM sequence (right end and A's when minus is
False, left end and T's when minus is True). A read antisense to the gene
has no run there and falls back to the alignment end.
"""

from tecap.constants import (
    NONTEMPLATED_MAX,
    SITE_RUN_MIN,
    SITE_RUN_TOLERANCE,
    TEMPLATED_MIN,
)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ALIGNED, _QUERY, _REF = (0, 7, 8), (1, 4), (2, 3)


def sense_3p_clip(read, minus):
    """Soft-clipped bases after the aligned 3' end, in transcript sense."""
    cig, seq = read.cigartuples, read.query_sequence
    if minus:
        n = cig[0][1] if cig[0][0] == 4 else 0
        return seq[:n].translate(_COMP)[::-1]
    n = cig[-1][1] if cig[-1][0] == 4 else 0
    return seq[len(seq) - n:] if n else ""


def inward_pairs(read, minus):
    """(query, ref) pairs of the aligned bases from the transcript 3' end
    inward, walked lazily from the CIGAR."""
    cig = read.cigartuples
    if minus:
        q, r, step = 0, read.reference_start, 1
    else:
        q, r, step, cig = read.query_length - 1, read.reference_end - 1, -1, reversed(cig)
    for op, n in cig:
        if op in _ALIGNED:
            for k in range(n):
                yield q + step * k, r + step * k
            q, r = q + step * n, r + step * n
        elif op in _QUERY:
            q += step * n
        elif op in _REF:
            r += step * n


def terminal_a_run(read, minus: bool, tolerance: int = SITE_RUN_TOLERANCE) -> tuple[int, int, int]:
    """Terminal A run around the aligned 3' end in transcript sense, with up to
    `tolerance` interior non-A bases.

    Grows rightwards over the 3' clip, then leftwards over the aligned bases,
    from one shared non-A budget; a non-A is taken only between A's. Returns
    (site, aligned_a, clip_a): the genomic boundary where the run starts (the
    alignment end when the run has no aligned base) and the A's in its aligned
    and clipped parts.
    """
    seq = read.query_sequence.translate(_COMP) if minus else read.query_sequence
    clip, pairs, left = sense_3p_clip(read, minus), inward_pairs(read, minus), []

    def aligned(t):
        """t-th aligned sense base from the 3' end as (base, ref), None past the 5' end."""
        while len(left) <= t:
            p = next(pairs, None)
            if p is None:
                return None
            left.append((seq[p[0]], p[1]))
        return left[t]

    budget, j, n = tolerance, 0, len(clip)
    while j < n:
        if clip[j] != "A":
            if not (budget and j + 1 < n and clip[j + 1] == "A" and (j or aligned(0)[0] == "A")):
                break
            budget -= 1
        j += 1
    t = aligned_a = 0
    while (cur := aligned(t)) is not None:
        if cur[0] == "A":
            aligned_a += 1
        else:
            nxt = aligned(t + 1)
            if not (budget and nxt is not None and nxt[0] == "A" and (t or j)):
                break
            budget -= 1
        t += 1
    if t == 0:
        site = read.reference_start if minus else read.reference_end
    else:
        site = left[t - 1][1] + 1 if minus else left[t - 1][1]
    return site, aligned_a, clip[:j].count("A")


def read_site(read, strand: str, anchor: str) -> tuple[int, bool]:
    """(site, observed) on the gene strand. alignment_end: the aligned 3' end,
    always observed. a_run_start: the start of the terminal A run when it has
    at least SITE_RUN_MIN A's, else the aligned 3' end and not observed. A
    record stored without SEQ has no run to read and falls back."""
    minus = strand == "-"
    end = read.reference_start if minus else read.reference_end
    if anchor == "alignment_end":
        return end, True
    if read.query_sequence is None:
        return end, False
    site, aligned_a, clip_a = terminal_a_run(read, minus)
    if aligned_a + clip_a >= SITE_RUN_MIN:
        return site, True
    return end, False


def label_of(a_count: int) -> str:
    """Label of the A count in the window downstream of a site."""
    if a_count >= TEMPLATED_MIN:
        return "templated"
    if a_count <= NONTEMPLATED_MAX:
        return "nontemplated"
    return "ambiguous"
