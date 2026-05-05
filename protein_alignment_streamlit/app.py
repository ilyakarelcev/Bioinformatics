from __future__ import annotations

import html
import re
import textwrap
from dataclasses import dataclass
from typing import Iterable

import streamlit as st
import streamlit.components.v1 as components

try:
    from Bio import BiopythonDeprecationWarning, pairwise2
    from Bio.Align import substitution_matrices
except Exception:  # pragma: no cover - handled in the UI
    pairwise2 = None
    substitution_matrices = None
    BiopythonDeprecationWarning = None

import warnings

if BiopythonDeprecationWarning is not None:
    warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


APP_TITLE = "Protein Pairwise Alignment Lab"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZ")
BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1
"""
AA_CLASS = {
    "A": "hydrophobic",
    "V": "hydrophobic",
    "I": "hydrophobic",
    "L": "hydrophobic",
    "M": "hydrophobic",
    "F": "aromatic",
    "Y": "aromatic",
    "W": "aromatic",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "C": "special",
    "G": "special",
    "P": "special",
    "K": "positive",
    "R": "positive",
    "H": "positive",
    "D": "negative",
    "E": "negative",
    "B": "ambiguous",
    "Z": "ambiguous",
    "X": "ambiguous",
}

EXAMPLE_FULL = (
    ">human_preproinsulin\n"
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)
EXAMPLE_FRAGMENT = (
    ">insulin_b_chain_like_fragment\n"
    "FVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGG"
)


@dataclass(frozen=True)
class AlignmentResult:
    algorithm: str
    seq1_name: str
    seq2_name: str
    raw_seq1: str
    raw_seq2: str
    aligned_seq1: str
    aligned_seq2: str
    score: float
    start: int
    end: int
    gap_open: float
    gap_extend: float


@dataclass(frozen=True)
class AlignmentMetrics:
    aligned_length: int
    paired_residues: int
    exact_matches: int
    positive_matches: int
    weak_matches: int
    mismatches: int
    gap_columns: int
    seq1_covered: int
    seq2_covered: int
    identity_alignment: float
    identity_paired: float
    similarity_alignment: float
    gap_percent: float
    seq1_coverage: float
    seq2_coverage: float


def parse_fasta_or_sequence(raw: str, fallback_name: str) -> tuple[str, str]:
    name = fallback_name
    body: list[str] = []
    for line in raw.splitlines():
        clean = line.strip()
        if not clean:
            continue
        if clean.startswith(">"):
            if name == fallback_name:
                name = clean[1:].strip() or fallback_name
            continue
        body.append(clean)
    sequence = re.sub(r"[^A-Za-z*.-]", "", "".join(body)).upper()
    return name, sequence


def normalize_protein_sequence(sequence: str) -> tuple[str, list[str]]:
    warnings_out: list[str] = []
    sequence = sequence.replace("-", "").replace(".", "").replace("*", "")
    replacements = {"U": "X", "O": "X", "J": "X"}
    normalized = []
    replaced: list[str] = []
    dropped: list[str] = []

    for char in sequence:
        if char in replacements:
            normalized.append(replacements[char])
            replaced.append(char)
        elif char in VALID_AA:
            normalized.append(char)
        else:
            dropped.append(char)

    if replaced:
        warnings_out.append(
            "Residues U/O/J are not represented in BLOSUM62 here and were converted to X: "
            + ", ".join(sorted(set(replaced)))
        )
    if dropped:
        warnings_out.append(
            "Unexpected symbols were removed: " + ", ".join(sorted(set(dropped)))
        )

    return "".join(normalized), warnings_out


@st.cache_resource(show_spinner=False)
def load_blosum62():
    if substitution_matrices is None:
        return parse_blosum62_text()
    return substitution_matrices.load("BLOSUM62")


@st.cache_resource(show_spinner=False)
def parse_blosum62_text() -> dict[tuple[str, str], int]:
    lines = [line.split() for line in BLOSUM62_TEXT.strip().splitlines()]
    header = lines[0]
    matrix: dict[tuple[str, str], int] = {}
    for row in lines[1:]:
        aa = row[0]
        for other, value in zip(header, row[1:]):
            matrix[(aa, other)] = int(value)
    return matrix


def matrix_score(matrix, aa1: str, aa2: str) -> float:
    if aa1 == "-" or aa2 == "-":
        return 0
    if isinstance(matrix, dict):
        return float(matrix.get((aa1, aa2), matrix.get((aa2, aa1), 0)))
    try:
        return float(matrix[aa1, aa2])
    except Exception:
        return 0


def fallback_pairwise_alignment(
    seq1: str,
    seq2: str,
    algorithm: str,
    gap_open: float,
    gap_extend: float,
) -> tuple[str, str, float, int, int]:
    matrix = load_blosum62()
    n = len(seq1)
    m = len(seq2)
    neg_inf = -10**12
    states = ("M", "X", "Y")

    score = {
        state: [[neg_inf] * (m + 1) for _ in range(n + 1)]
        for state in states
    }
    pointer: dict[str, list[list[tuple[str, int, int] | None]]] = {
        state: [[None] * (m + 1) for _ in range(n + 1)]
        for state in states
    }

    score["M"][0][0] = 0
    if algorithm == "Needleman-Wunsch":
        for i in range(1, n + 1):
            score["X"][i][0] = gap_open + (i - 1) * gap_extend
            pointer["X"][i][0] = ("X", i - 1, 0) if i > 1 else ("M", 0, 0)
        for j in range(1, m + 1):
            score["Y"][0][j] = gap_open + (j - 1) * gap_extend
            pointer["Y"][0][j] = ("Y", 0, j - 1) if j > 1 else ("M", 0, 0)
    else:
        for state in states:
            for i in range(n + 1):
                score[state][i][0] = 0
            for j in range(m + 1):
                score[state][0][j] = 0

    best_state = "M"
    best_i = n
    best_j = m
    best_score = neg_inf

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            subst = matrix_score(matrix, seq1[i - 1], seq2[j - 1])

            candidates_m = [
                (score["M"][i - 1][j - 1] + subst, "M"),
                (score["X"][i - 1][j - 1] + subst, "X"),
                (score["Y"][i - 1][j - 1] + subst, "Y"),
            ]
            best_m, state_m = max(candidates_m, key=lambda item: item[0])
            score["M"][i][j] = best_m
            pointer["M"][i][j] = (state_m, i - 1, j - 1)

            candidates_x = [
                (score["M"][i - 1][j] + gap_open, "M"),
                (score["X"][i - 1][j] + gap_extend, "X"),
                (score["Y"][i - 1][j] + gap_open, "Y"),
            ]
            best_x, state_x = max(candidates_x, key=lambda item: item[0])
            score["X"][i][j] = best_x
            pointer["X"][i][j] = (state_x, i - 1, j)

            candidates_y = [
                (score["M"][i][j - 1] + gap_open, "M"),
                (score["X"][i][j - 1] + gap_open, "X"),
                (score["Y"][i][j - 1] + gap_extend, "Y"),
            ]
            best_y, state_y = max(candidates_y, key=lambda item: item[0])
            score["Y"][i][j] = best_y
            pointer["Y"][i][j] = (state_y, i, j - 1)

            if algorithm == "Smith-Waterman":
                for state in states:
                    if score[state][i][j] < 0:
                        score[state][i][j] = 0
                        pointer[state][i][j] = None
                    if score[state][i][j] > best_score:
                        best_score = score[state][i][j]
                        best_state = state
                        best_i = i
                        best_j = j

    if algorithm == "Needleman-Wunsch":
        terminal = [(score[state][n][m], state) for state in states]
        best_score, best_state = max(terminal, key=lambda item: item[0])
        best_i = n
        best_j = m

    aligned_1: list[str] = []
    aligned_2: list[str] = []
    state = best_state
    i = best_i
    j = best_j
    end_column = 0

    while i > 0 or j > 0:
        if algorithm == "Smith-Waterman" and score[state][i][j] <= 0:
            break
        previous = pointer[state][i][j]
        if previous is None:
            break

        prev_state, prev_i, prev_j = previous
        if prev_i == i - 1 and prev_j == j - 1:
            aligned_1.append(seq1[i - 1])
            aligned_2.append(seq2[j - 1])
        elif prev_i == i - 1 and prev_j == j:
            aligned_1.append(seq1[i - 1])
            aligned_2.append("-")
        else:
            aligned_1.append("-")
            aligned_2.append(seq2[j - 1])

        state = prev_state
        i = prev_i
        j = prev_j
        end_column += 1

    aligned_1.reverse()
    aligned_2.reverse()
    return "".join(aligned_1), "".join(aligned_2), float(best_score), 0, end_column


def run_alignment(
    seq1_name: str,
    seq2_name: str,
    seq1: str,
    seq2: str,
    algorithm: str,
    gap_open: float,
    gap_extend: float,
) -> AlignmentResult:
    matrix = load_blosum62()
    if matrix is None:
        raise RuntimeError("BLOSUM62 matrix could not be loaded.")

    if pairwise2 is None:
        aligned_seq1, aligned_seq2, score, start, end = fallback_pairwise_alignment(
            seq1,
            seq2,
            algorithm,
            gap_open,
            gap_extend,
        )
    else:
        aligner = pairwise2.align.localds if algorithm == "Smith-Waterman" else pairwise2.align.globalds
        alignments = aligner(
            seq1,
            seq2,
            matrix,
            gap_open,
            gap_extend,
            one_alignment_only=True,
        )
        if not alignments:
            raise RuntimeError("No alignment was produced for these sequences.")

        best = alignments[0]
        aligned_seq1 = best.seqA
        aligned_seq2 = best.seqB
        score = float(best.score)
        start = int(best.start)
        end = int(best.end)

        if algorithm == "Smith-Waterman":
            aligned_seq1 = aligned_seq1[start:end]
            aligned_seq2 = aligned_seq2[start:end]

    return AlignmentResult(
        algorithm=algorithm,
        seq1_name=seq1_name,
        seq2_name=seq2_name,
        raw_seq1=seq1,
        raw_seq2=seq2,
        aligned_seq1=aligned_seq1,
        aligned_seq2=aligned_seq2,
        score=score,
        start=start,
        end=end,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def analyze_alignment(result: AlignmentResult) -> AlignmentMetrics:
    matrix = load_blosum62()
    aligned_length = len(result.aligned_seq1)
    paired = exact = positive = weak = mismatches = gaps = 0
    seq1_covered = seq2_covered = 0

    for aa1, aa2 in zip(result.aligned_seq1, result.aligned_seq2):
        if aa1 != "-":
            seq1_covered += 1
        if aa2 != "-":
            seq2_covered += 1
        if aa1 == "-" or aa2 == "-":
            gaps += 1
            continue

        paired += 1
        score = matrix_score(matrix, aa1, aa2)
        if aa1 == aa2:
            exact += 1
            positive += 1
        elif score > 0:
            positive += 1
        elif score == 0:
            weak += 1
        else:
            mismatches += 1

    safe_aligned = max(aligned_length, 1)
    safe_paired = max(paired, 1)
    return AlignmentMetrics(
        aligned_length=aligned_length,
        paired_residues=paired,
        exact_matches=exact,
        positive_matches=positive,
        weak_matches=weak,
        mismatches=mismatches,
        gap_columns=gaps,
        seq1_covered=seq1_covered,
        seq2_covered=seq2_covered,
        identity_alignment=exact / safe_aligned * 100,
        identity_paired=exact / safe_paired * 100,
        similarity_alignment=positive / safe_aligned * 100,
        gap_percent=gaps / safe_aligned * 100,
        seq1_coverage=seq1_covered / max(len(result.raw_seq1), 1) * 100,
        seq2_coverage=seq2_covered / max(len(result.raw_seq2), 1) * 100,
    )


def connector_symbol(matrix, aa1: str, aa2: str) -> str:
    if aa1 == "-" or aa2 == "-":
        return ""
    if aa1 == aa2:
        return "|"
    score = matrix_score(matrix, aa1, aa2)
    if score > 0:
        return ":"
    if score == 0:
        return "."
    return ""


def residue_class(aa: str, other: str, color_mode: str, matrix) -> str:
    if aa == "-":
        return "gap"
    if color_mode == "By residue property":
        return AA_CLASS.get(aa, "ambiguous")
    if other == "-":
        return "gap-neighbor"
    if aa == other:
        return "exact"
    score = matrix_score(matrix, aa, other)
    if score > 0:
        return "positive-match"
    if score == 0:
        return "weak-match"
    return "mismatch"


def position_labels(aligned: str) -> list[str]:
    labels: list[str] = []
    position = 0
    for aa in aligned:
        if aa == "-":
            labels.append("")
            continue
        position += 1
        labels.append(str(position) if position == 1 or position % 10 == 0 else "")
    return labels


def chunk_ranges(length: int, size: int) -> Iterable[range]:
    for start in range(0, length, size):
        yield range(start, min(start + size, length))


def build_alignment_html(result: AlignmentResult, color_mode: str, wrap: int) -> str:
    matrix = load_blosum62()
    labels_1 = position_labels(result.aligned_seq1)
    labels_2 = position_labels(result.aligned_seq2)
    chunks: list[str] = []

    for indexes in chunk_ranges(len(result.aligned_seq1), wrap):
        cols = len(indexes)
        row_top: list[str] = []
        row_a: list[str] = []
        row_mid: list[str] = []
        row_b: list[str] = []
        row_bottom: list[str] = []

        for idx in indexes:
            aa1 = result.aligned_seq1[idx]
            aa2 = result.aligned_seq2[idx]
            symbol = connector_symbol(matrix, aa1, aa2)
            class_1 = residue_class(aa1, aa2, color_mode, matrix)
            class_2 = residue_class(aa2, aa1, color_mode, matrix)
            title_1 = html.escape(f"{result.seq1_name}: {aa1} vs {aa2}")
            title_2 = html.escape(f"{result.seq2_name}: {aa2} vs {aa1}")
            row_top.append(f'<div class="pos">{html.escape(labels_1[idx])}</div>')
            row_a.append(f'<div class="tile {class_1}" title="{title_1}">{html.escape(aa1)}</div>')
            row_mid.append(f'<div class="connector">{html.escape(symbol)}</div>')
            row_b.append(f'<div class="tile {class_2}" title="{title_2}">{html.escape(aa2)}</div>')
            row_bottom.append(f'<div class="pos">{html.escape(labels_2[idx])}</div>')

        grid_style = f"--cols:{cols};"
        chunks.append(
            f"""
            <section class="alignment-chunk">
              <div class="chunk-label">
                <span>{html.escape(result.seq1_name)}</span>
                <span>columns {indexes.start + 1}-{indexes.stop}</span>
                <span>{html.escape(result.seq2_name)}</span>
              </div>
              <div class="grid" style="{grid_style}">{''.join(row_top)}</div>
              <div class="grid" style="{grid_style}">{''.join(row_a)}</div>
              <div class="grid connectors" style="{grid_style}">{''.join(row_mid)}</div>
              <div class="grid" style="{grid_style}">{''.join(row_b)}</div>
              <div class="grid" style="{grid_style}">{''.join(row_bottom)}</div>
            </section>
            """
        )

    return f"""
    <style>
      :root {{
        color-scheme: light;
      }}
      body {{
        margin: 0;
        font-family: Inter, Segoe UI, Arial, sans-serif;
        background: #f7f8fb;
        color: #172033;
      }}
      .viewer {{
        padding: 18px;
        border: 1px solid #d9dee8;
        background: #ffffff;
        border-radius: 8px;
      }}
      .legend {{
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        margin-bottom: 16px;
        font-size: 12px;
      }}
      .legend-item {{
        display: inline-flex;
        align-items: center;
        gap: 6px;
        padding: 5px 8px;
        border: 1px solid #d9dee8;
        border-radius: 999px;
        background: #fff;
      }}
      .swatch {{
        width: 14px;
        height: 14px;
        border-radius: 4px;
        display: inline-block;
        border: 1px solid rgba(23, 32, 51, 0.12);
      }}
      .alignment-chunk {{
        overflow-x: auto;
        padding: 12px 0 16px;
        border-top: 1px solid #eef1f6;
      }}
      .alignment-chunk:first-of-type {{
        border-top: 0;
      }}
      .chunk-label {{
        min-width: max-content;
        display: grid;
        grid-template-columns: 1fr auto 1fr;
        gap: 14px;
        margin-bottom: 7px;
        color: #687289;
        font-size: 12px;
        font-weight: 650;
      }}
      .chunk-label span:last-child {{
        text-align: right;
      }}
      .grid {{
        display: grid;
        grid-template-columns: repeat(var(--cols), 26px);
        gap: 3px;
        min-width: max-content;
      }}
      .tile {{
        height: 26px;
        display: flex;
        align-items: center;
        justify-content: center;
        border-radius: 6px;
        font-family: Consolas, Menlo, monospace;
        font-weight: 800;
        font-size: 14px;
        border: 1px solid rgba(23, 32, 51, 0.09);
        box-sizing: border-box;
      }}
      .connector {{
        height: 16px;
        display: flex;
        align-items: center;
        justify-content: center;
        font-family: Consolas, Menlo, monospace;
        font-weight: 900;
        color: #526078;
      }}
      .connectors {{
        margin: 1px 0;
      }}
      .pos {{
        height: 14px;
        display: flex;
        align-items: center;
        justify-content: center;
        color: #8a94a8;
        font-size: 9px;
        font-family: Consolas, Menlo, monospace;
      }}
      .exact {{ background: #10b981; color: #06281d; }}
      .positive-match {{ background: #8bd3ff; color: #093047; }}
      .weak-match {{ background: #ffe08a; color: #493700; }}
      .mismatch {{ background: #e7eaf0; color: #394357; }}
      .gap, .gap-neighbor {{ background: #ffd1d1; color: #7a1f1f; }}
      .hydrophobic {{ background: #d8df8f; color: #2f360b; }}
      .polar {{ background: #a7e7c0; color: #123822; }}
      .positive {{ background: #8bbcff; color: #0d2d5c; }}
      .negative {{ background: #ff9f9f; color: #5d1111; }}
      .aromatic {{ background: #cdb7ff; color: #2f1b69; }}
      .special {{ background: #f4c178; color: #4b2c00; }}
      .ambiguous {{ background: #d7dce5; color: #354052; }}
    </style>
    <div class="viewer">
      <div class="legend">
        <span class="legend-item"><span class="swatch exact"></span> exact match</span>
        <span class="legend-item"><span class="swatch positive-match"></span> positive BLOSUM62</span>
        <span class="legend-item"><span class="swatch weak-match"></span> neutral BLOSUM62</span>
        <span class="legend-item"><span class="swatch mismatch"></span> mismatch</span>
        <span class="legend-item"><span class="swatch gap"></span> gap</span>
      </div>
      {''.join(chunks)}
    </div>
    """


def interpret(result: AlignmentResult, metrics: AlignmentMetrics) -> list[str]:
    notes: list[str] = []
    if result.algorithm == "Smith-Waterman":
        notes.append(
            "Local alignment is useful here because it focuses on the best matching region, "
            "so a peptide, domain, or partial protein can still align cleanly."
        )
        if min(metrics.seq1_coverage, metrics.seq2_coverage) < 55 <= max(metrics.seq1_coverage, metrics.seq2_coverage):
            notes.append(
                "The coverage is asymmetric: one sequence may be a fragment of the other or represent a shared domain."
            )
    else:
        notes.append(
            "Global alignment compares both sequences end to end, which is best when they are expected to be homologous across most of their length."
        )

    if metrics.identity_paired >= 90:
        notes.append("Very high identity: the proteins are probably very close variants or nearly the same sequence.")
    elif metrics.identity_paired >= 50:
        notes.append("High identity across paired residues: this strongly suggests close relationship between the aligned regions.")
    elif metrics.identity_paired >= 25 and metrics.similarity_alignment >= 35:
        notes.append("Moderate identity with additional conservative substitutions: this can be compatible with more distant homology.")
    else:
        notes.append("Low identity: the alignment alone gives weak evidence for relatedness, especially if coverage is also low.")

    if metrics.gap_percent >= 20:
        notes.append("Many gap columns are present, which may indicate insertions, deletions, or different protein architecture.")
    elif metrics.gap_percent <= 5:
        notes.append("Few gaps are present, so the aligned region has a fairly consistent length.")

    notes.append("These statements are alignment-based clues, not a proof of homology or function.")
    return notes


def make_report(result: AlignmentResult, metrics: AlignmentMetrics, notes: list[str]) -> str:
    connector = "".join(
        connector_symbol(load_blosum62(), aa1, aa2)
        for aa1, aa2 in zip(result.aligned_seq1, result.aligned_seq2)
    )
    wrapped_blocks: list[str] = []
    width = 80
    for start in range(0, len(result.aligned_seq1), width):
        stop = start + width
        wrapped_blocks.append(result.aligned_seq1[start:stop])
        wrapped_blocks.append(connector[start:stop])
        wrapped_blocks.append(result.aligned_seq2[start:stop])
        wrapped_blocks.append("")

    return "\n".join(
        [
            APP_TITLE,
            "",
            f"Algorithm: {result.algorithm}",
            "Matrix: BLOSUM62",
            f"Gap open: {result.gap_open}",
            f"Gap extension: {result.gap_extend}",
            f"Score: {result.score:.2f}",
            "",
            f"Sequence 1: {result.seq1_name}, length {len(result.raw_seq1)} aa",
            f"Sequence 2: {result.seq2_name}, length {len(result.raw_seq2)} aa",
            f"Alignment length: {metrics.aligned_length}",
            f"Identity over alignment: {metrics.identity_alignment:.1f}%",
            f"Identity over paired residues: {metrics.identity_paired:.1f}%",
            f"Similarity over alignment: {metrics.similarity_alignment:.1f}%",
            f"Gap columns: {metrics.gap_columns} ({metrics.gap_percent:.1f}%)",
            f"Coverage sequence 1: {metrics.seq1_coverage:.1f}%",
            f"Coverage sequence 2: {metrics.seq2_coverage:.1f}%",
            "",
            "Interpretation:",
            *[f"- {note}" for note in notes],
            "",
            "Alignment:",
            *wrapped_blocks,
        ]
    )


def inject_page_css() -> None:
    st.markdown(
        """
        <style>
          .stApp {
            background:
              linear-gradient(180deg, #f6f8fb 0%, #ffffff 34%);
          }
          [data-testid="stMetric"] {
            background: #ffffff;
            border: 1px solid #dde3ee;
            border-radius: 8px;
            padding: 12px 14px;
          }
          [data-testid="stMetricLabel"] {
            color: #5d6880;
          }
          div[data-testid="stForm"] {
            border: 1px solid #dde3ee;
            border-radius: 8px;
            background: #ffffff;
          }
        </style>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    st.set_page_config(page_title=APP_TITLE, page_icon="AA", layout="wide")
    inject_page_css()

    st.title(APP_TITLE)
    st.caption(
        "Compare two amino-acid sequences with BLOSUM62, Smith-Waterman or Needleman-Wunsch, "
        "then inspect the alignment as a readable biological map."
    )

    if pairwise2 is None:
        st.info(
            "Biopython is not installed in this environment, so the app is using its built-in "
            "BLOSUM62 alignment engine. Install the requirements to use Biopython's implementation."
        )

    with st.sidebar:
        st.header("Alignment settings")
        algorithm_label = st.radio(
            "Algorithm",
            [
                "Smith-Waterman: local",
                "Needleman-Wunsch: global",
            ],
            index=0,
        )
        algorithm = "Smith-Waterman" if algorithm_label.startswith("Smith") else "Needleman-Wunsch"
        gap_open = st.number_input("Gap open penalty", value=-10.0, min_value=-50.0, max_value=0.0, step=0.5)
        gap_extend = st.number_input("Gap extension penalty", value=-0.5, min_value=-20.0, max_value=0.0, step=0.1)
        color_mode = st.radio(
            "Tile colors",
            ["By alignment result", "By residue property"],
            index=0,
            horizontal=True,
        )
        wrap = st.slider("Columns per block", min_value=30, max_value=120, value=70, step=10)

        st.divider()
        st.markdown(
            textwrap.dedent(
                """
                **Default scoring**

                Matrix: `BLOSUM62`

                Exact match: shown as `|`

                Positive substitution: shown as `:`

                Neutral substitution: shown as `.`
                """
            )
        )

    with st.form("alignment_form"):
        col1, col2 = st.columns(2)
        with col1:
            raw_1 = st.text_area("Sequence 1", value=EXAMPLE_FULL, height=220)
        with col2:
            raw_2 = st.text_area("Sequence 2", value=EXAMPLE_FRAGMENT, height=220)

        submitted = st.form_submit_button("Align sequences", type="primary", use_container_width=True)

    if not submitted:
        st.info("Paste FASTA or plain amino-acid sequences, choose settings, then run the alignment.")
        st.stop()

    seq1_name, raw_seq1 = parse_fasta_or_sequence(raw_1, "Sequence 1")
    seq2_name, raw_seq2 = parse_fasta_or_sequence(raw_2, "Sequence 2")
    seq1, warnings_1 = normalize_protein_sequence(raw_seq1)
    seq2, warnings_2 = normalize_protein_sequence(raw_seq2)

    for message in warnings_1 + warnings_2:
        st.warning(message)

    if not seq1 or not seq2:
        st.error("Both sequences must contain at least one amino-acid residue.")
        st.stop()

    with st.spinner("Aligning sequences..."):
        result = run_alignment(seq1_name, seq2_name, seq1, seq2, algorithm, gap_open, gap_extend)
        metrics = analyze_alignment(result)
        notes = interpret(result, metrics)

    metric_cols = st.columns(6)
    metric_cols[0].metric("Score", f"{result.score:.1f}")
    metric_cols[1].metric("Identity", f"{metrics.identity_alignment:.1f}%")
    metric_cols[2].metric("Paired identity", f"{metrics.identity_paired:.1f}%")
    metric_cols[3].metric("Similarity", f"{metrics.similarity_alignment:.1f}%")
    metric_cols[4].metric("Gaps", f"{metrics.gap_percent:.1f}%")
    metric_cols[5].metric("Length", metrics.aligned_length)

    tab_alignment, tab_metrics, tab_interpretation, tab_export = st.tabs(
        ["Alignment map", "Metrics", "Interpretation", "Export"]
    )

    with tab_alignment:
        viewer_height = min(900, max(360, (metrics.aligned_length // max(wrap, 1) + 1) * 150))
        components.html(
            build_alignment_html(result, color_mode or "By alignment result", wrap),
            height=viewer_height,
            scrolling=True,
        )

    with tab_metrics:
        left, right = st.columns(2)
        with left:
            st.subheader("Alignment")
            st.write(f"Algorithm: **{result.algorithm}**")
            st.write("Substitution matrix: **BLOSUM62**")
            st.write(f"Gap open / extension: **{gap_open} / {gap_extend}**")
            st.write(f"Aligned columns: **{metrics.aligned_length}**")
            st.write(f"Paired residue columns: **{metrics.paired_residues}**")
            st.write(f"Exact matches: **{metrics.exact_matches}**")
            st.write(f"Positive BLOSUM62 columns: **{metrics.positive_matches}**")
            st.write(f"Neutral columns: **{metrics.weak_matches}**")
            st.write(f"Mismatches: **{metrics.mismatches}**")
            st.write(f"Gap columns: **{metrics.gap_columns}**")
        with right:
            st.subheader("Coverage")
            st.write(f"{seq1_name}: **{metrics.seq1_covered} / {len(seq1)} aa** ({metrics.seq1_coverage:.1f}%)")
            st.write(f"{seq2_name}: **{metrics.seq2_covered} / {len(seq2)} aa** ({metrics.seq2_coverage:.1f}%)")
            if result.algorithm == "Smith-Waterman":
                st.write(f"Local alignment span in gapped alignment: **{result.start + 1}-{result.end}**")
            st.write(f"Raw sequence 1 length: **{len(seq1)} aa**")
            st.write(f"Raw sequence 2 length: **{len(seq2)} aa**")

    with tab_interpretation:
        st.subheader("What the alignment suggests")
        for note in notes:
            st.markdown(f"- {note}")

    with tab_export:
        report = make_report(result, metrics, notes)
        st.download_button(
            "Download TXT report",
            data=report,
            file_name="protein_pairwise_alignment_report.txt",
            mime="text/plain",
            use_container_width=True,
        )
        st.code(report, language="text")


if __name__ == "__main__":
    main()
