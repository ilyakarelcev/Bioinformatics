from __future__ import annotations

import html
import json
import re
import textwrap
import warnings
from dataclasses import dataclass
from typing import Iterable

import streamlit as st
import streamlit.components.v1 as components

warnings.filterwarnings("ignore", message="Bio.pairwise2 has been deprecated*")

try:
    from Bio import BiopythonDeprecationWarning, pairwise2
    from Bio.Align import substitution_matrices
except Exception:  # pragma: no cover - handled in the UI
    pairwise2 = None
    substitution_matrices = None
    BiopythonDeprecationWarning = None

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
AA_NAMES = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamic acid",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "B": "Asparagine or aspartic acid",
    "Z": "Glutamine or glutamic acid",
    "X": "Unknown amino acid",
    "-": "Gap",
}

EXAMPLE_FULL = (
    ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n"
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)
EXAMPLE_FRAGMENT = (
    ">sp|P01315|INS_PIG Insulin OS=Sus scrofa OX=9823 GN=INS PE=1 SV=2\n"
    "MALWTRLLPLLALLALWAPAPAQAFVNQHLCGSHLVEALYLVCGERGFFYTPKARREAEN"
    "PQAGAVELGGGLGGLQALALEGPPQKRGIVEQCCTSICSLYQLENYCN"
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


def connector_explanation(matrix, aa1: str, aa2: str) -> tuple[str, str]:
    if aa1 == "-" or aa2 == "-":
        return "Gap", "One sequence has an insertion/deletion here."
    if aa1 == aa2:
        return "Exact match", "`|` means the same amino acid in both sequences."
    score = matrix_score(matrix, aa1, aa2)
    if score > 0:
        return "Conservative substitution", "`:` means BLOSUM62 scores this replacement as favorable."
    if score == 0:
        return "Weak similarity", "`.` means neutral/weak similarity in BLOSUM62."
    return "Mismatch", "No symbol means an unfavorable substitution."


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


def build_overview_html(result: AlignmentResult) -> str:
    matrix = load_blosum62()
    cells: list[str] = []
    for idx, (aa1, aa2) in enumerate(zip(result.aligned_seq1, result.aligned_seq2), start=1):
        if aa1 == "-" or aa2 == "-":
            class_name = "gap"
        elif aa1 == aa2:
            class_name = "exact"
        elif matrix_score(matrix, aa1, aa2) > 0:
            class_name = "positive-match"
        elif matrix_score(matrix, aa1, aa2) == 0:
            class_name = "weak-match"
        else:
            class_name = "mismatch"
        cells.append(f'<span class="overview-cell {class_name}" title="Column {idx}: {aa1}/{aa2}"></span>')

    return f"""
    <style>
      body {{
        margin: 0;
        font-family: Inter, Segoe UI, Arial, sans-serif;
        background: transparent;
      }}
      .overview {{
        border: 1px solid #d9dee8;
        background: #ffffff;
        border-radius: 8px;
        padding: 12px;
      }}
      .overview-labels {{
        display: flex;
        justify-content: space-between;
        color: #687289;
        font-size: 12px;
        font-weight: 700;
        margin-bottom: 8px;
      }}
      .overview-track {{
        display: grid;
        grid-template-columns: repeat({len(cells)}, minmax(2px, 1fr));
        gap: 1px;
        min-width: 100%;
        height: 34px;
        overflow: hidden;
        border-radius: 6px;
        border: 1px solid #edf0f5;
        background: #f5f7fb;
      }}
      .overview-cell {{
        min-width: 2px;
        height: 100%;
      }}
      .exact {{ background: #10b981; }}
      .positive-match {{ background: #8bd3ff; }}
      .weak-match {{ background: #ffe08a; }}
      .mismatch {{ background: #d4dae5; }}
      .gap {{ background: #ff9f9f; }}
    </style>
    <div class="overview">
      <div class="overview-labels">
        <span>1</span>
        <span>Alignment overview: {len(cells)} columns</span>
        <span>{len(cells)}</span>
      </div>
      <div class="overview-track">{''.join(cells)}</div>
    </div>
    """


def build_alignment_html(
    result: AlignmentResult,
    color_mode: str,
    wrap: int,
    start_column: int = 1,
    end_column: int | None = None,
) -> str:
    matrix = load_blosum62()
    labels_1 = position_labels(result.aligned_seq1)
    labels_2 = position_labels(result.aligned_seq2)
    chunks: list[str] = []
    alignment_length = len(result.aligned_seq1)
    start_index = max(0, start_column - 1)
    end_index = alignment_length if end_column is None else min(alignment_length, end_column)
    selected_length = max(0, end_index - start_index)

    for local_indexes in chunk_ranges(selected_length, wrap):
        indexes = range(start_index + local_indexes.start, start_index + local_indexes.stop)
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
        overflow-x: auto;
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


def build_interactive_alignment_html(result: AlignmentResult, color_mode: str, app_theme: str) -> str:
    matrix = load_blosum62()
    seq1_position = 0
    seq2_position = 0
    columns: list[dict[str, str | int]] = []

    for idx, (aa1, aa2) in enumerate(zip(result.aligned_seq1, result.aligned_seq2), start=1):
        if aa1 != "-":
            seq1_position += 1
        if aa2 != "-":
            seq2_position += 1

        symbol = connector_symbol(matrix, aa1, aa2)
        relation_title, relation_detail = connector_explanation(matrix, aa1, aa2)
        columns.append(
            {
                "column": idx,
                "aa1": aa1,
                "aa2": aa2,
                "aa1Name": AA_NAMES.get(aa1, "Unknown residue"),
                "aa2Name": AA_NAMES.get(aa2, "Unknown residue"),
                "relationTitle": relation_title,
                "relationDetail": relation_detail,
                "seq1Pos": seq1_position if aa1 != "-" else "",
                "seq2Pos": seq2_position if aa2 != "-" else "",
                "symbol": symbol,
                "overviewClass": residue_class(aa1, aa2, "By alignment result", matrix),
                "class1": residue_class(aa1, aa2, color_mode, matrix),
                "class2": residue_class(aa2, aa1, color_mode, matrix),
            }
        )

    payload = json.dumps(
        {
            "seq1Name": result.seq1_name,
            "seq2Name": result.seq2_name,
            "columns": columns,
        }
    )
    is_dark = app_theme == "Dark"
    viewer_colors = {
        "page_bg": "#080814" if is_dark else "#ffffff",
        "panel_bg": "#111122" if is_dark else "#ffffff",
        "panel_border": "#2c2a4a" if is_dark else "#d0d7de",
        "text": "#f4f2ff" if is_dark else "#24292f",
        "muted": "#a8a2cf" if is_dark else "#57606a",
        "soft_border": "#282747" if is_dark else "#d8dee4",
        "scale": "#8f89c6" if is_dark else "#57606a",
        "tick": "#d7d2ff" if is_dark else "#24292f",
        "track_bg": "#0b0b19" if is_dark else "#f6f8fa",
        "track_border": "#302d55" if is_dark else "#d0d7de",
        "viewport_bg": "rgba(124, 92, 255, 0.20)" if is_dark else "rgba(31, 136, 61, 0.14)",
        "viewport_border": "rgba(164, 143, 255, 0.88)" if is_dark else "rgba(31, 136, 61, 0.66)",
        "handle_bg": "rgba(20, 18, 42, 0.94)" if is_dark else "rgba(255, 255, 255, 0.9)",
        "position": "#a8a2cf" if is_dark else "#57606a",
        "connector": "#d7d2ff" if is_dark else "#57606a",
        "hover": "#7c5cff" if is_dark else "#0969da",
        "paired": "#2dd4bf" if is_dark else "#9a6700",
    }
    tile_colors = {
        "exact_bg": "#7c5cff" if is_dark else "#34d058",
        "exact_fg": "#f8f6ff" if is_dark else "#052e16",
        "positive_bg": "#2dd4bf" if is_dark else "#8be9df",
        "positive_fg": "#041f20" if is_dark else "#063b3b",
        "weak_bg": "#f0abfc" if is_dark else "#ffd33d",
        "weak_fg": "#321044" if is_dark else "#3b2300",
        "mismatch_bg": "#26243d" if is_dark else "#eaeef2",
        "mismatch_fg": "#d8d3ff" if is_dark else "#24292f",
        "gap_bg": "#fb7185" if is_dark else "#ffebe9",
        "gap_fg": "#2b0710" if is_dark else "#82071e",
        "hydrophobic_bg": "#575174" if is_dark else "#d8dee4",
        "hydrophobic_fg": "#f2efff" if is_dark else "#24292f",
        "polar_bg": "#38bdf8" if is_dark else "#a5d6a7",
        "polar_fg": "#061b2b" if is_dark else "#143d1a",
        "positive_charge_bg": "#818cf8" if is_dark else "#79c0ff",
        "positive_charge_fg": "#101342" if is_dark else "#05264c",
        "negative_bg": "#f472b6" if is_dark else "#ffaba8",
        "negative_fg": "#331023" if is_dark else "#5f0f0b",
        "aromatic_bg": "#c084fc" if is_dark else "#d2a8ff",
        "aromatic_fg": "#25103f" if is_dark else "#3b0f70",
        "special_bg": "#fbbf24" if is_dark else "#ffdf8b",
        "special_fg": "#2e1b00" if is_dark else "#442900",
        "ambiguous_bg": "#4b4868" if is_dark else "#d0d7de",
        "ambiguous_fg": "#ebe7ff" if is_dark else "#24292f",
    }

    return f"""
    <style>
      body {{
        margin: 0;
        font-family: Inter, Segoe UI, Arial, sans-serif;
        background: {viewer_colors["page_bg"]};
        color: {viewer_colors["text"]};
      }}
      .alignment-browser {{
        border: 1px solid {viewer_colors["panel_border"]};
        border-radius: 8px;
        padding: 14px;
        background: {viewer_colors["panel_bg"]};
        box-sizing: border-box;
        user-select: none;
      }}
      .browser-head {{
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 12px;
        margin-bottom: 10px;
        color: {viewer_colors["muted"]};
        font-size: 12px;
        font-weight: 700;
      }}
      .overview {{
        position: relative;
        height: 68px;
        padding: 0 0 14px;
      }}
      .scale-line {{
        position: absolute;
        left: 0;
        right: 0;
        top: 10px;
        height: 1px;
        background: {viewer_colors["scale"]};
      }}
      .tick {{
        position: absolute;
        top: 10px;
        width: 1px;
        height: 10px;
        background: {viewer_colors["scale"]};
      }}
      .tick span {{
        position: absolute;
        top: 12px;
        left: 50%;
        transform: translateX(-50%);
        font-size: 10px;
        color: {viewer_colors["tick"]};
        white-space: nowrap;
      }}
      .track {{
        position: absolute;
        left: 0;
        right: 0;
        bottom: 0;
        height: 28px;
        display: grid;
        grid-template-columns: repeat(var(--cols), minmax(1px, 1fr));
        overflow: hidden;
        background: {viewer_colors["track_bg"]};
        border: 1px solid {viewer_colors["track_border"]};
        border-radius: 3px;
      }}
      .overview-cell {{
        min-width: 1px;
        height: 100%;
      }}
      .viewport {{
        position: absolute;
        top: 10px;
        bottom: 0;
        background: {viewer_colors["viewport_bg"]};
        border: 1px solid {viewer_colors["viewport_border"]};
        cursor: grab;
        box-sizing: border-box;
        z-index: 4;
      }}
      .viewport:active {{
        cursor: grabbing;
      }}
      .handle {{
        position: absolute;
        top: -2px;
        bottom: -2px;
        width: 8px;
        background: {viewer_colors["handle_bg"]};
        border: 1px solid {viewer_colors["viewport_border"]};
        box-sizing: border-box;
        cursor: ew-resize;
      }}
      .handle.left {{
        left: -5px;
      }}
      .handle.right {{
        right: -5px;
      }}
      .range-readout {{
        margin: 10px 0 8px;
        color: {viewer_colors["muted"]};
        font-size: 12px;
        font-weight: 700;
      }}
      .detail {{
        border-top: 1px solid {viewer_colors["soft_border"]};
        padding-top: 10px;
      }}
      .lane {{
        position: relative;
        overflow: hidden;
        width: 100%;
        box-sizing: border-box;
      }}
      .lane-inner {{
        display: grid;
        grid-template-columns: repeat(var(--cols), var(--tile-width));
        gap: var(--tile-gap);
        width: max-content;
        transform: translateX(var(--offset-x));
        will-change: transform;
      }}
      .detail-scale {{
        height: 16px;
        margin-bottom: 3px;
      }}
      .detail-pos {{
        text-align: center;
        color: {viewer_colors["position"]};
        font-size: 9px;
        font-family: Consolas, Menlo, monospace;
        overflow: hidden;
        white-space: nowrap;
      }}
      .tile {{
        position: relative;
        height: 28px;
        width: var(--tile-width);
        min-width: var(--tile-width);
        display: flex;
        align-items: center;
        justify-content: center;
        border-radius: 4px;
        border: 1px solid rgba(23, 32, 51, 0.09);
        box-sizing: border-box;
        font-family: Consolas, Menlo, monospace;
        font-size: var(--tile-font-size);
        font-weight: 800;
        overflow: hidden;
        white-space: nowrap;
      }}
      .tile.is-hovered {{
        outline: 2px solid #f8fafc;
        outline-offset: -2px;
        box-shadow: 0 0 0 2px {viewer_colors["hover"]}, 0 0 18px rgba(88, 166, 255, 0.35);
        z-index: 3;
      }}
      .tile.is-paired {{
        outline: 2px solid {viewer_colors["paired"]};
        outline-offset: -2px;
        box-shadow: 0 0 0 2px rgba(210, 153, 34, 0.55);
        z-index: 2;
      }}
      .connector-row {{
        height: 22px;
        align-items: center;
      }}
      .connector {{
        min-width: 0;
        height: 22px;
        display: flex;
        align-items: center;
        justify-content: center;
        overflow: hidden;
      }}
      .connector.is-paired {{
        background: rgba(250, 204, 21, 0.12);
      }}
      .connector-mark {{
        display: block;
        flex: 0 0 auto;
      }}
      .connector-mark.exact-line {{
        width: 2px;
        height: 18px;
        border-radius: 999px;
        background: {viewer_colors["connector"]};
      }}
      .connector-mark.positive-dots {{
        position: relative;
        width: 6px;
        height: 18px;
      }}
      .connector-mark.positive-dots::before,
      .connector-mark.positive-dots::after {{
        content: "";
        position: absolute;
        left: 50%;
        width: 4px;
        height: 4px;
        border-radius: 999px;
        background: {viewer_colors["connector"]};
        transform: translateX(-50%);
      }}
      .connector-mark.positive-dots::before {{
        top: 4px;
      }}
      .connector-mark.positive-dots::after {{
        bottom: 4px;
      }}
      .connector-mark.weak-dot {{
        width: 4px;
        height: 4px;
        border-radius: 999px;
        background: {viewer_colors["connector"]};
      }}
      .compact .tile span {{
        transform: scale(0.78);
      }}
      .hide-residue-letters .tile span {{
        opacity: 0;
      }}
      .micro .tile {{
        height: 22px;
        border-radius: 0;
        border-left: 0;
        border-right: 0;
      }}
      .name-row {{
        display: grid;
        grid-template-columns: 120px 1fr;
        gap: 8px;
        align-items: center;
        margin: 3px 0;
      }}
      .seq-name {{
        color: {viewer_colors["muted"]};
        font-size: 12px;
        font-weight: 800;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }}
      .tooltip {{
        position: fixed;
        z-index: 30;
        display: none;
        max-width: 260px;
        padding: 9px 10px;
        border-radius: 7px;
        background: rgba(17, 24, 39, 0.96);
        color: #ffffff;
        box-shadow: 0 8px 24px rgba(15, 23, 42, 0.22);
        font-size: 12px;
        line-height: 1.35;
        pointer-events: none;
      }}
      .tooltip strong {{
        display: block;
        margin-bottom: 2px;
        font-size: 13px;
      }}
      .tooltip span {{
        color: #cbd5e1;
      }}
      .exact {{ background: {tile_colors["exact_bg"]}; color: {tile_colors["exact_fg"]}; }}
      .positive-match {{ background: {tile_colors["positive_bg"]}; color: {tile_colors["positive_fg"]}; }}
      .weak-match {{ background: {tile_colors["weak_bg"]}; color: {tile_colors["weak_fg"]}; }}
      .mismatch {{ background: {tile_colors["mismatch_bg"]}; color: {tile_colors["mismatch_fg"]}; }}
      .gap, .gap-neighbor {{ background: {tile_colors["gap_bg"]}; color: {tile_colors["gap_fg"]}; }}
      .hydrophobic {{ background: {tile_colors["hydrophobic_bg"]}; color: {tile_colors["hydrophobic_fg"]}; }}
      .polar {{ background: {tile_colors["polar_bg"]}; color: {tile_colors["polar_fg"]}; }}
      .positive {{ background: {tile_colors["positive_charge_bg"]}; color: {tile_colors["positive_charge_fg"]}; }}
      .negative {{ background: {tile_colors["negative_bg"]}; color: {tile_colors["negative_fg"]}; }}
      .aromatic {{ background: {tile_colors["aromatic_bg"]}; color: {tile_colors["aromatic_fg"]}; }}
      .special {{ background: {tile_colors["special_bg"]}; color: {tile_colors["special_fg"]}; }}
      .ambiguous {{ background: {tile_colors["ambiguous_bg"]}; color: {tile_colors["ambiguous_fg"]}; }}
    </style>

    <div class="alignment-browser" id="browser">
      <div class="browser-head">
        <span id="leftLabel"></span>
        <span id="rightLabel"></span>
      </div>
      <div class="overview" id="overview">
        <div class="scale-line"></div>
        <div id="ticks"></div>
        <div class="track" id="track"></div>
        <div class="viewport" id="viewport">
          <span class="handle left" data-mode="left"></span>
          <span class="handle right" data-mode="right"></span>
        </div>
      </div>
      <div class="range-readout" id="readout"></div>
      <div class="detail" id="detail">
        <div class="name-row">
          <div class="seq-name"></div>
          <div class="lane detail-scale"><div class="lane-inner" id="detailScale"></div></div>
        </div>
        <div class="name-row">
          <div class="seq-name" id="seq1Name"></div>
          <div class="lane tile-row"><div class="lane-inner" id="row1"></div></div>
        </div>
        <div class="name-row">
          <div class="seq-name"></div>
          <div class="lane connector-row"><div class="lane-inner" id="connectors"></div></div>
        </div>
        <div class="name-row">
          <div class="seq-name" id="seq2Name"></div>
          <div class="lane tile-row"><div class="lane-inner" id="row2"></div></div>
        </div>
      </div>
      <div class="tooltip" id="tooltip"></div>
    </div>

    <script>
      const payload = {payload};
      const columns = payload.columns;
      const total = columns.length;
      const browser = document.getElementById("browser");
      const overview = document.getElementById("overview");
      const track = document.getElementById("track");
      const viewport = document.getElementById("viewport");
      const ticks = document.getElementById("ticks");
      const readout = document.getElementById("readout");
      const detail = document.getElementById("detail");
      const detailScale = document.getElementById("detailScale");
      const row1 = document.getElementById("row1");
      const row2 = document.getElementById("row2");
      const connectors = document.getElementById("connectors");
      const seq1Name = document.getElementById("seq1Name");
      const seq2Name = document.getElementById("seq2Name");
      const leftLabel = document.getElementById("leftLabel");
      const rightLabel = document.getElementById("rightLabel");
      const tooltip = document.getElementById("tooltip");

      seq1Name.textContent = payload.seq1Name;
      seq2Name.textContent = payload.seq2Name;
      leftLabel.textContent = "1";
      rightLabel.textContent = String(total);

      const initialSize = total <= 160 ? total : Math.min(120, total);
      let start = 0.0;
      let end = initialSize;
      let drag = null;
      const minCols = Math.min(8, total);

      function escapeHtml(value) {{
        return String(value)
          .replaceAll("&", "&amp;")
          .replaceAll("<", "&lt;")
          .replaceAll(">", "&gt;")
          .replaceAll('"', "&quot;");
      }}

      function chooseTickStep(n) {{
        if (n <= 120) return 10;
        if (n <= 400) return 25;
        if (n <= 1000) return 50;
        if (n <= 3000) return 100;
        return 500;
      }}

      function columnFromClientX(clientX) {{
        const rect = overview.getBoundingClientRect();
        const ratio = Math.min(1, Math.max(0, (clientX - rect.left) / rect.width));
        return ratio * total;
      }}

      function clampRange() {{
        start = Math.max(0, Math.min(start, total - minCols));
        end = Math.max(start + minCols, Math.min(end, total));
      }}

      function updateViewport() {{
        clampRange();
        viewport.style.left = (start / total * 100) + "%";
        viewport.style.width = ((end - start) / total * 100) + "%";
        readout.textContent = `Columns ${{Math.floor(start) + 1}}-${{Math.ceil(end)}} of ${{total}}`;
      }}

      function tileHtml(col, aaKey, classKey, posKey, otherKey) {{
        const aa = col[aaKey];
        const other = col[otherKey];
        const nameKey = aaKey === "aa1" ? "aa1Name" : "aa2Name";
        const otherNameKey = aaKey === "aa1" ? "aa2Name" : "aa1Name";
        return `<div class="tile ${{col[classKey]}}"
          data-index="${{col.column - 1}}"
          data-aa="${{escapeHtml(aa)}}"
          data-name="${{escapeHtml(col[nameKey])}}"
          data-other-name="${{escapeHtml(col[otherNameKey])}}"
          data-symbol="${{escapeHtml(col.symbol || "none")}}"
          data-relation-title="${{escapeHtml(col.relationTitle)}}"
          data-relation-detail="${{escapeHtml(col.relationDetail)}}"
          data-position="${{escapeHtml(col[posKey] || "gap")}}"
          data-column="${{col.column}}"
          data-other="${{escapeHtml(other)}}"><span>${{escapeHtml(aa)}}</span></div>`;
      }}

      function connectorHtml(col) {{
        let mark = "";
        if (col.symbol === "|") mark = '<span class="connector-mark exact-line"></span>';
        else if (col.symbol === ":") mark = '<span class="connector-mark positive-dots"></span>';
        else if (col.symbol === ".") mark = '<span class="connector-mark weak-dot"></span>';
        return `<div class="connector" data-index="${{col.column - 1}}">${{mark}}</div>`;
      }}

      function renderDetail() {{
        const visibleSpan = Math.max(minCols, end - start);
        const laneWidth = Math.max(1, row1.parentElement.getBoundingClientRect().width);
        let roughTileWidth = laneWidth / visibleSpan;
        let gap = roughTileWidth < 8 ? 0 : 2;
        let tileWidth = (laneWidth - gap * Math.max(0, visibleSpan - 1)) / visibleSpan;
        if (tileWidth < 1) {{
          gap = 0;
          tileWidth = laneWidth / visibleSpan;
        }}
        const offset = -(start * (tileWidth + gap));
        detail.style.setProperty("--cols", total);
        detail.style.setProperty("--tile-width", tileWidth + "px");
        detail.style.setProperty("--tile-gap", gap + "px");
        detail.style.setProperty("--offset-x", offset + "px");
        detail.style.setProperty("--tile-font-size", Math.max(7, Math.min(13, tileWidth * 0.82)) + "px");
        detail.style.setProperty("--connector-font-size", Math.max(7, Math.min(12, tileWidth * 0.9)) + "px");
        detail.classList.toggle("compact", tileWidth < 10);
        detail.classList.toggle("hide-residue-letters", tileWidth < 6);
        detail.classList.toggle("micro", tileWidth < 5);

        detailScale.innerHTML = columns.map((col) => {{
          const step = tileWidth >= 18 ? 5 : tileWidth >= 9 ? 10 : tileWidth >= 4 ? 20 : 50;
          const show = col.column === 1 || col.column === total || col.column % step === 0 || visibleSpan <= 35;
          return `<div class="detail-pos">${{show ? col.column : ""}}</div>`;
        }}).join("");
        row1.innerHTML = columns.map((col) => tileHtml(col, "aa1", "class1", "seq1Pos", "aa2")).join("");
        connectors.innerHTML = columns.map((col) => connectorHtml(col)).join("");
        row2.innerHTML = columns.map((col) => tileHtml(col, "aa2", "class2", "seq2Pos", "aa1")).join("");
      }}

      function renderOverview() {{
        track.style.setProperty("--cols", total);
        track.innerHTML = columns.map((col) => `<span class="overview-cell ${{col.overviewClass}}"></span>`).join("");
        const step = chooseTickStep(total);
        const tickValues = [1];
        for (let value = step; value < total; value += step) tickValues.push(value);
        if (!tickValues.includes(total)) tickValues.push(total);
        ticks.innerHTML = tickValues.map((value) => {{
          const left = total === 1 ? 0 : ((value - 1) / (total - 1) * 100);
          return `<div class="tick" style="left:${{left}}%"><span>${{value}}</span></div>`;
        }}).join("");
      }}

      function rerender() {{
        updateViewport();
        renderDetail();
      }}

      viewport.addEventListener("pointerdown", (event) => {{
        event.preventDefault();
        viewport.setPointerCapture(event.pointerId);
        const mode = event.target.dataset.mode || "move";
        drag = {{
          mode,
          x: event.clientX,
          start,
          end,
        }};
      }});

      viewport.addEventListener("pointermove", (event) => {{
        if (!drag) return;
        const deltaCols = columnFromClientX(event.clientX) - columnFromClientX(drag.x);
        if (drag.mode === "left") {{
          start = Math.min(drag.end - minCols, Math.max(0, drag.start + deltaCols));
          end = drag.end;
        }} else if (drag.mode === "right") {{
          start = drag.start;
          end = Math.max(drag.start + minCols, Math.min(total, drag.end + deltaCols));
        }} else {{
          const width = drag.end - drag.start;
          start = Math.max(0, Math.min(total - width, drag.start + deltaCols));
          end = start + width;
        }}
        rerender();
      }});

      viewport.addEventListener("pointerup", () => {{
        drag = null;
      }});

      overview.addEventListener("pointerdown", (event) => {{
        if (event.target === viewport || event.target.classList.contains("handle")) return;
        const width = end - start;
        const center = columnFromClientX(event.clientX);
        start = Math.max(0, Math.min(total - width, center - width / 2));
        end = start + width;
        rerender();
      }});

      window.addEventListener("resize", renderDetail);
      function clearHighlight() {{
        detail.querySelectorAll(".is-hovered, .is-paired").forEach((node) => {{
          node.classList.remove("is-hovered", "is-paired");
        }});
      }}

      function highlightColumn(tile) {{
        clearHighlight();
        const index = Number(tile.dataset.index);
        const first = row1.children[index];
        const second = row2.children[index];
        const link = connectors.children[index];
        tile.classList.add("is-hovered");
        if (first && first !== tile) first.classList.add("is-paired");
        if (second && second !== tile) second.classList.add("is-paired");
        if (link) link.classList.add("is-paired");
      }}

      detail.addEventListener("pointerover", (event) => {{
        const tile = event.target.closest(".tile");
        if (!tile) return;
        highlightColumn(tile);
        tooltip.innerHTML = `<strong>${{tile.dataset.aa}} - ${{tile.dataset.name}}</strong>
          <span>${{tile.dataset.sequence}}</span><br>
          Position: ${{tile.dataset.position}} · Alignment column: ${{tile.dataset.column}}<br>
          Compared with: ${{tile.dataset.other}}`;
        tooltip.innerHTML = `<strong>${{tile.dataset.aa}} - ${{tile.dataset.name}}</strong>
          <span>Position ${{tile.dataset.position}} · column ${{tile.dataset.column}}</span><br>
          Compared with: ${{tile.dataset.other}} - ${{tile.dataset.otherName}}<br>
          Symbol: ${{tile.dataset.symbol}} · ${{tile.dataset.relationTitle}}<br>
          <span>${{tile.dataset.relationDetail}}</span>`;
        tooltip.style.display = "block";
      }});
      detail.addEventListener("pointermove", (event) => {{
        if (tooltip.style.display !== "block") return;
        const padding = 12;
        const rect = tooltip.getBoundingClientRect();
        let left = event.clientX - rect.width / 2;
        let top = event.clientY - rect.height - 14;
        left = Math.max(padding, Math.min(window.innerWidth - rect.width - padding, left));
        if (top < padding) top = event.clientY + 16;
        tooltip.style.left = left + "px";
        tooltip.style.top = top + "px";
      }});
      detail.addEventListener("pointerleave", () => {{
        clearHighlight();
        tooltip.style.display = "none";
      }});
      renderOverview();
      rerender();
    </script>
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


def inject_page_css(app_theme: str) -> None:
    is_dark = app_theme == "Dark"
    colors = {
        "app_bg": (
            "radial-gradient(circle at 18% -8%, rgba(124, 92, 255, 0.22) 0%, rgba(8, 8, 20, 0) 38%), linear-gradient(180deg, #080814 0%, #05050d 100%)"
            if is_dark
            else "linear-gradient(180deg, #f6f8fa 0%, #ffffff 42%)"
        ),
        "panel_bg": "#111122" if is_dark else "#ffffff",
        "panel_bg_soft": "#0b0b19" if is_dark else "#f6f8fa",
        "panel_border": "#2c2a4a" if is_dark else "#d0d7de",
        "metric_label": "#a8a2cf" if is_dark else "#57606a",
        "metric_value": "#f4f2ff" if is_dark else "#24292f",
        "text": "#f4f2ff" if is_dark else "#24292f",
        "muted": "#a8a2cf" if is_dark else "#57606a",
        "textarea_bg": "#0b0b19" if is_dark else "#ffffff",
        "textarea_border": "#343056" if is_dark else "#8c959f",
        "textarea_focus": "#8b5cf6" if is_dark else "#1f883d",
        "button_bg": "#7c5cff" if is_dark else "#1f883d",
        "button_hover": "#8b7cff" if is_dark else "#1a7f37",
        "button_text": "#ffffff",
        "tab_active": "#a78bfa" if is_dark else "#1f883d",
        "code_bg": "#0b0b19" if is_dark else "#f6f8fa",
    }
    st.markdown(
        f"""
        <style>
          .stApp {{
            background: {colors["app_bg"]};
            color: {colors["text"]};
          }}
          .stApp, .stApp p, .stApp span, .stApp label {{
            color: {colors["text"]};
          }}
          .stCaptionContainer, .stCaptionContainer p, small {{
            color: {colors["muted"]} !important;
          }}
          [data-testid="stMetric"] {{
            background: {colors["panel_bg"]};
            border: 1px solid {colors["panel_border"]};
            border-radius: 8px;
            padding: 12px 14px;
          }}
          [data-testid="stMetricLabel"] {{
            color: {colors["metric_label"]};
          }}
          [data-testid="stMetricValue"] {{
            color: {colors["metric_value"]};
          }}
          [data-testid="stMetricDelta"] {{
            color: {colors["muted"]};
          }}
          [data-testid="stWidgetLabel"],
          [data-testid="stWidgetLabel"] p,
          [data-testid="stTextArea"] label,
          [data-testid="stTextArea"] label p {{
            color: {colors["text"]} !important;
            font-weight: 650;
          }}
          textarea {{
            background: {colors["textarea_bg"]} !important;
            border-color: {colors["textarea_border"]} !important;
            color: {colors["text"]} !important;
            caret-color: {colors["textarea_focus"]} !important;
          }}
          textarea:focus {{
            border-color: {colors["textarea_focus"]} !important;
            box-shadow: 0 0 0 1px {colors["textarea_focus"]} !important;
          }}
          input, div[data-baseweb="input"] input {{
            color: {colors["text"]} !important;
          }}
          div[data-baseweb="input"] {{
            background: {colors["textarea_bg"]} !important;
            border-color: {colors["textarea_border"]} !important;
          }}
          .stButton > button,
          .stDownloadButton > button {{
            border-radius: 8px;
            border: 1px solid {colors["panel_border"]};
            color: {colors["text"]};
            background: {colors["panel_bg"]};
          }}
          .stButton > button[kind="primary"],
          .stButton > button[data-testid="baseButton-primary"] {{
            background: {colors["button_bg"]} !important;
            border-color: {colors["button_bg"]} !important;
            color: {colors["button_text"]} !important;
          }}
          .stButton > button[kind="primary"]:hover,
          .stButton > button[data-testid="baseButton-primary"]:hover {{
            background: {colors["button_hover"]} !important;
            border-color: {colors["button_hover"]} !important;
            color: {colors["button_text"]} !important;
          }}
          .stButton > button:disabled {{
            opacity: 0.55;
            color: {colors["muted"]} !important;
            background: {colors["panel_bg_soft"]} !important;
          }}
          div[data-testid="stExpander"], div[data-testid="stAlert"] {{
            border-color: {colors["panel_border"]};
          }}
          div[data-testid="stTabs"] [role="tablist"] {{
            border-bottom-color: {colors["panel_border"]};
          }}
          div[data-testid="stTabs"] button[role="tab"] p {{
            color: {colors["muted"]} !important;
          }}
          div[data-testid="stTabs"] button[role="tab"][aria-selected="true"] p {{
            color: {colors["tab_active"]} !important;
          }}
          div[data-testid="stTabs"] [data-baseweb="tab-highlight"] {{
            background-color: {colors["tab_active"]} !important;
          }}
          section[data-testid="stSidebar"] {{
            border-right: 1px solid {colors["panel_border"]};
            background: {colors["panel_bg"]};
          }}
          div[data-testid="stForm"] {{
            border: 1px solid {colors["panel_border"]};
            border-radius: 8px;
            background: {colors["panel_bg"]};
          }}
          code, pre {{
            background: {colors["code_bg"]} !important;
            color: {colors["text"]} !important;
          }}
          div[data-testid="stMarkdownContainer"] {{
            color: {colors["text"]};
          }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    st.set_page_config(page_title=APP_TITLE, page_icon="AA", layout="wide")

    with st.sidebar:
        app_theme = st.radio("Theme", ["Dark", "Light"], index=0, horizontal=True)
        st.divider()

    inject_page_css(app_theme)

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

    col1, col2 = st.columns(2)
    with col1:
        raw_1 = st.text_area("Sequence 1", value=EXAMPLE_FULL, height=220, key="sequence_1_input")
    with col2:
        raw_2 = st.text_area("Sequence 2", value=EXAMPLE_FRAGMENT, height=220, key="sequence_2_input")

    run_col, auto_col = st.columns([3, 1])
    with auto_col:
        auto_align = st.checkbox(
            "Auto Align",
            value=False,
            help="Recalculate automatically whenever sequence text or settings change.",
        )
    with run_col:
        submitted = st.button(
            "Align sequences",
            type="primary",
            use_container_width=True,
            disabled=auto_align,
        )
    if auto_align:
        st.caption("Auto Align is on: sequence and scoring changes recalculate immediately.")

    should_align = auto_align or submitted

    if should_align:
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

        st.session_state["alignment_payload"] = {
            "result": result,
            "metrics": metrics,
            "notes": notes,
            "seq1_name": seq1_name,
            "seq2_name": seq2_name,
            "seq1": seq1,
            "seq2": seq2,
            "gap_open": gap_open,
            "gap_extend": gap_extend,
        }
    elif "alignment_payload" not in st.session_state:
        st.info("Paste FASTA or plain amino-acid sequences, choose settings, then run the alignment.")
        st.stop()

    payload = st.session_state["alignment_payload"]
    result = payload["result"]
    metrics = payload["metrics"]
    notes = payload["notes"]
    seq1_name = payload["seq1_name"]
    seq2_name = payload["seq2_name"]
    seq1 = payload["seq1"]
    seq2 = payload["seq2"]
    used_gap_open = payload["gap_open"]
    used_gap_extend = payload["gap_extend"]

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
        components.html(
            build_interactive_alignment_html(result, color_mode or "By alignment result", app_theme),
            height=330,
            scrolling=False,
        )

    with tab_metrics:
        left, right = st.columns(2)
        with left:
            st.subheader("Alignment")
            st.write(f"Algorithm: **{result.algorithm}**")
            st.write("Substitution matrix: **BLOSUM62**")
            st.write(f"Gap open / extension: **{used_gap_open} / {used_gap_extend}**")
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
