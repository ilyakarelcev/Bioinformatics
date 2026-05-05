from __future__ import annotations

import json
import re
from typing import Any

import pandas as pd
import requests
import streamlit as st


BASE_URL = "https://rest.uniprot.org"
DOCS = "https://www.uniprot.org/api-documentation"

STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {codon for codon, aa in STANDARD_CODE.items() if aa == "*"}
START_CODONS = {"ATG"}

DEFAULT_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "length",
    "mass",
    "reviewed",
    "cc_function",
    "xref_pdb",
    "xref_alphafold",
]

EXAMPLE_PROTEINS = [
    {
        "id": "P01308",
        "name": "Инсулин человека",
        "note": "маленький гормон, но очень подробно аннотирован",
    },
    {
        "id": "P00698",
        "name": "Лизоцим курицы",
        "note": "классический небольшой фермент, удобен для простого JSON",
    },
    {
        "id": "P69905",
        "name": "Гемоглобин альфа человека",
        "note": "короткая субъединица, типичная Swiss-Prot запись",
    },
    {
        "id": "P68871",
        "name": "Гемоглобин бета человека",
        "note": "похож на предыдущий, но с другой аннотацией и вариантами",
    },
    {
        "id": "P0A7V8",
        "name": "Рибосомный белок S4 E. coli",
        "note": "бактериальный белок, хороший не-человеческий пример",
    },
    {
        "id": "P0A6F5",
        "name": "Chaperonin GroEL E. coli",
        "note": "покрупнее и с понятной функцией шаперона",
    },
    {
        "id": "Q8N158",
        "name": "Glypican-2 человека",
        "note": "человеческий белок, обычно менее на слуху, чем инсулин",
    },
    {
        "id": "A0A0B4J2F0",
        "name": "PIOS1 человека",
        "note": "пример TrEMBL/unreviewed записи: данных обычно меньше",
    },
]

VIEW_OPTIONS = {
    "Карточка": "card",
    "JSON-дерево": "json_tree",
    "JSON простым текстом": "json_plain",
    "FASTA": "fasta",
    "Компактная таблица TSV": "tsv",
}


st.set_page_config(
    page_title="Лаборатория UniProt API",
    page_icon="UP",
    layout="wide",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Outfit', sans-serif !important;
}

/* App background gradient */
.stApp {
    background: radial-gradient(circle at top left, #1f1816 0%, #0d0b0a 100%) !important;
}

/* Stylish headers */
h1, h2, h3 {
    color: #f05b32 !important;
    text-shadow: 0 0 15px rgba(240, 91, 50, 0.2);
    font-weight: 700 !important;
    letter-spacing: -0.5px;
}

/* Metric styling */
[data-testid="stMetricValue"] {
    color: #ff8c42 !important;
    font-weight: 700 !important;
    text-shadow: 0 0 10px rgba(255, 140, 66, 0.3);
}
[data-testid="stMetricLabel"] {
    color: #c7b8b2 !important;
    font-weight: 600 !important;
}

/* Primary buttons */
button[data-testid="baseButton-primary"] {
    background: linear-gradient(135deg, #f05b32 0%, #d43b11 100%) !important;
    border: none !important;
    box-shadow: 0 4px 15px rgba(240, 91, 50, 0.3) !important;
    transition: all 0.3s ease !important;
    border-radius: 8px !important;
    color: #ffffff !important;
    font-weight: 600 !important;
}
button[data-testid="baseButton-primary"]:hover {
    box-shadow: 0 6px 20px rgba(240, 91, 50, 0.5) !important;
    transform: translateY(-2px);
}

/* Secondary buttons */
button[data-testid="baseButton-secondary"] {
    border: 1px solid rgba(240, 91, 50, 0.4) !important;
    color: #ff8c42 !important;
    transition: all 0.3s ease !important;
    background: rgba(26, 22, 20, 0.6) !important;
    border-radius: 8px !important;
    font-weight: 600 !important;
}
button[data-testid="baseButton-secondary"]:hover {
    border: 1px solid #f05b32 !important;
    box-shadow: 0 0 15px rgba(240, 91, 50, 0.2) !important;
    color: #f05b32 !important;
}

/* Borders and containers */
[data-testid="stVerticalBlockBorderWrapper"] {
    border: 1px solid rgba(255, 140, 66, 0.15) !important;
    border-radius: 16px !important;
    background: rgba(41, 36, 33, 0.4) !important;
    backdrop-filter: blur(12px) !important;
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4) !important;
}

/* Input fields */
.stTextInput>div>div>input, .stTextArea>div>div>textarea {
    background-color: rgba(13, 11, 10, 0.5) !important;
    border: 1px solid rgba(240, 91, 50, 0.3) !important;
    color: #f4f0ec !important;
    border-radius: 8px !important;
    transition: all 0.3s ease !important;
}
.stTextInput>div>div>input:focus, .stTextArea>div>div>textarea:focus {
    border: 1px solid #f05b32 !important;
    box-shadow: 0 0 10px rgba(240, 91, 50, 0.2) !important;
}

/* Tabs */
[data-testid="stTabs"] button {
    font-weight: 600 !important;
    font-size: 1.05rem !important;
    color: #c7b8b2 !important;
}
[data-testid="stTabs"] button[aria-selected="true"] {
    color: #f05b32 !important;
}
[data-testid="stTabs"] button[aria-selected="true"] > div[data-testid="stMarkdownContainer"] > p {
    color: #f05b32 !important;
}

/* Code blocks */
.stCodeBlock {
    border: 1px solid rgba(255, 140, 66, 0.15) !important;
    border-radius: 8px !important;
    background: #110e0c !important;
}

/* Selectbox */
[data-testid="stSelectbox"] > div > div {
    background-color: rgba(13, 11, 10, 0.5) !important;
    border: 1px solid rgba(240, 91, 50, 0.3) !important;
    border-radius: 8px !important;
}
[data-testid="stSelectbox"] > div > div:focus-within {
    border: 1px solid #f05b32 !important;
    box-shadow: 0 0 10px rgba(240, 91, 50, 0.2) !important;
}

/* Multiselect */
[data-testid="stMultiSelect"] > div > div {
    background-color: rgba(13, 11, 10, 0.5) !important;
    border: 1px solid rgba(240, 91, 50, 0.3) !important;
    border-radius: 8px !important;
}
[data-testid="stMultiSelect"] > div > div:focus-within {
    border: 1px solid #f05b32 !important;
    box-shadow: 0 0 10px rgba(240, 91, 50, 0.2) !important;
}

/* Sliders */
.stSlider div[data-testid="stTickBar"] > div {
    background-color: #f05b32 !important;
}

/* Expander headers */
.streamlit-expanderHeader {
    font-weight: 600 !important;
    color: #f05b32 !important;
}
</style>
""", unsafe_allow_html=True)



@st.cache_data(ttl=60 * 30, show_spinner=False)
def get_text(url: str, params: dict[str, str] | None = None) -> str:
    response = requests.get(url, params=params, timeout=25)
    response.raise_for_status()
    return response.text


@st.cache_data(ttl=60 * 30, show_spinner=False)
def get_json(url: str, params: dict[str, str] | None = None) -> dict[str, Any]:
    response = requests.get(url, params=params, timeout=25)
    response.raise_for_status()
    return response.json()


def normalize_accession(value: str) -> str:
    return value.strip().split()[0].upper()


def clean_nucleotide_sequence(raw: str) -> str:
    lines = [line.strip() for line in raw.splitlines() if not line.strip().startswith(">")]
    return re.sub(r"[^ACGTUacgtu]", "", "".join(lines)).upper().replace("U", "T")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate_dna(sequence: str) -> str:
    amino_acids = []
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        amino_acids.append(STANDARD_CODE.get(codon, "X"))
    return "".join(amino_acids)


def remove_nested_orfs(orfs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    kept: list[dict[str, Any]] = []
    for orf in sorted(orfs, key=lambda item: item["nt_length"], reverse=True):
        nested = any(
            orf["strand"] == other["strand"]
            and orf["start"] >= other["start"]
            and orf["end"] <= other["end"]
            for other in kept
        )
        if not nested:
            kept.append(orf)
    return sorted(kept, key=lambda item: item["nt_length"], reverse=True)


def find_orfs(
    sequence: str,
    min_nt_length: int = 90,
    start_mode: str = "ATG only",
    include_partial: bool = False,
    ignore_nested: bool = True,
) -> list[dict[str, Any]]:
    sequence = clean_nucleotide_sequence(sequence)
    if len(sequence) < 3:
        return []

    orfs: list[dict[str, Any]] = []
    strands = [("+", sequence), ("-", reverse_complement(sequence))]
    seq_len = len(sequence)

    for strand, scan_sequence in strands:
        for frame in range(3):
            index = frame
            while index <= len(scan_sequence) - 3:
                codon = scan_sequence[index:index + 3]
                is_start = codon in START_CODONS if start_mode == "ATG only" else codon not in STOP_CODONS
                if not is_start:
                    index += 3
                    continue

                stop_index = None
                for cursor in range(index + 3, len(scan_sequence) - 2, 3):
                    if scan_sequence[cursor:cursor + 3] in STOP_CODONS:
                        stop_index = cursor
                        break

                if stop_index is None and not include_partial:
                    index += 3
                    continue

                end_index = stop_index + 3 if stop_index is not None else len(scan_sequence) - ((len(scan_sequence) - index) % 3)
                nt_sequence = scan_sequence[index:end_index]
                if len(nt_sequence) < min_nt_length:
                    index += 3
                    continue

                if strand == "+":
                    start = index + 1
                    end = end_index
                    frame_label = f"+{frame + 1}"
                else:
                    start = seq_len - end_index + 1
                    end = seq_len - index
                    frame_label = f"-{frame + 1}"

                protein = translate_dna(nt_sequence)
                if protein.endswith("*"):
                    protein = protein[:-1]

                orfs.append(
                    {
                        "frame": frame_label,
                        "strand": strand,
                        "start": start,
                        "end": end,
                        "nt_length": len(nt_sequence),
                        "aa_length": len(protein),
                        "start_codon": codon,
                        "stop_codon": scan_sequence[stop_index:stop_index + 3] if stop_index is not None else "partial",
                        "protein": protein,
                        "nucleotide": nt_sequence,
                    }
                )
                index += 3 if start_mode == "ATG only" else len(nt_sequence)

    if ignore_nested:
        orfs = remove_nested_orfs(orfs)

    return sorted(orfs, key=lambda item: item["nt_length"], reverse=True)


def parse_ids_with_scores(raw: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line in raw.splitlines():
        clean = line.strip()
        if not clean or clean.startswith("#"):
            continue

        parts = re.split(r"[\s,;]+", clean)
        accession = normalize_accession(parts[0])
        score = parts[1] if len(parts) > 1 else ""
        if accession:
            rows.append({"accession": accession, "score": score})
    return rows


def record_summary(record: dict[str, Any]) -> dict[str, Any]:
    genes = record.get("genes") or []
    gene_names = []
    for gene in genes:
        primary = gene.get("geneName", {}).get("value")
        if primary:
            gene_names.append(primary)

    sequence = record.get("sequence") or {}
    protein_description = record.get("proteinDescription") or {}
    recommended = protein_description.get("recommendedName") or {}
    full_name = recommended.get("fullName", {}).get("value")

    if not full_name:
        submission_names = protein_description.get("submissionNames") or []
        if submission_names:
            full_name = submission_names[0].get("fullName", {}).get("value")

    organism = record.get("organism") or {}
    comments = record.get("comments") or []
    function_texts = []
    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            for text in comment.get("texts") or []:
                value = text.get("value")
                if value:
                    function_texts.append(value)

    return {
        "accession": record.get("primaryAccession", ""),
        "entry": record.get("uniProtkbId", ""),
        "protein": full_name or "Нет рекомендуемого названия белка",
        "genes": ", ".join(gene_names) or "Нет названия гена",
        "organism": organism.get("scientificName", "Неизвестный организм"),
        "taxon": organism.get("taxonId", ""),
        "length": sequence.get("length", ""),
        "mass": sequence.get("molWeight", ""),
        "reviewed": bool(record.get("entryType", "").lower().startswith("uniprotkb reviewed")),
        "sequence": sequence.get("value", ""),
        "function": "\n\n".join(function_texts),
    }


def show_summary_card(summary: dict[str, Any], score: str | None = None) -> None:
    reviewed = "Проверенная запись / Swiss-Prot" if summary["reviewed"] else "Непроверенная запись / TrEMBL"
    score_label = f"Сходство: {score}" if score else "Сходство не передано"

    with st.container(border=True):
        left, right = st.columns([2, 1])
        with left:
            st.subheader(f"{summary['accession']} · {summary['protein']}")
            st.caption(f"{summary['entry']} · {reviewed}")
            st.write(summary["function"] or "В этом ответе нет комментария о функции.")
        with right:
            st.metric("Длина", f"{summary['length']} aa" if summary["length"] else "н/д")
            st.metric("Масса", f"{summary['mass']} Da" if summary["mass"] else "н/д")
            st.write(f"**Ген:** {summary['genes']}")
            st.write(f"**Организм:** {summary['organism']}")
            st.write(f"**Taxon ID:** {summary['taxon'] or 'н/д'}")
            st.write(f"**{score_label}**")


def fetch_record(accession: str) -> dict[str, Any]:
    return get_json(f"{BASE_URL}/uniprotkb/{accession}.json")


def fetch_fasta(accession: str) -> str:
    return get_text(f"{BASE_URL}/uniprotkb/{accession}.fasta")


def fetch_flat_fields(accession: str, fields: list[str]) -> pd.DataFrame:
    text = get_text(
        f"{BASE_URL}/uniprotkb/search",
        params={
            "query": f"accession:{accession}",
            "format": "tsv",
            "fields": ",".join(fields),
            "size": "1",
        },
    )
    rows = [line.split("\t") for line in text.strip().splitlines()]
    if len(rows) < 2:
        return pd.DataFrame()
    return pd.DataFrame(rows[1:], columns=rows[0])


def api_error(error: Exception) -> None:
    st.error("Запрос к UniProt API не удался.")
    st.exception(error)


st.title("Лаборатория UniProt API")
st.caption("Небольшой Streamlit-интерфейс для просмотра белков, последовательностей и карточек похожих кандидатов из UniProt.")

with st.expander("Зачем нужен этот прототип", expanded=True):
    st.markdown(
        """
        Большая система, похоже, будет работать примерно так:

        1. Пользователь вводит аминокислотную последовательность и вопрос обычным текстом.
        2. Бэкенд отделяет последовательность от вопроса.
        3. Поиск по графу или эмбеддингам находит похожие белки.
        4. Этот поиск возвращает один или несколько UniProt ID, часто еще и со score похожести.
        5. Этот интерфейс обращается к UniProt и превращает ID в понятные данные: название белка, ген, организм, функцию, последовательность, ссылки на структуры и другие базы.

        В этой схеме UniProt не ищет похожесть сам. Он скорее слой метаданных и доказательств вокруг ID, которые уже нашла ваша система.
        """
    )

tab_lookup, tab_search, tab_hits, tab_sequence, tab_orf, tab_docs = st.tabs(
    [
        "Поиск по ID",
        "Поиск в UniProtKB",
        "Карточки похожих",
        "Ввод последовательности",
        "ORF-анализ ДНК",
        "Карта API",
    ]
)

with tab_lookup:
    st.header("Получить запись по UniProt ID")
    st.write("Эта вкладка нужна, когда UniProt accession уже известен, например его вернул поиск по эмбеддингам.")

    example_labels = [
        f"{item['id']} — {item['name']}: {item['note']}" for item in EXAMPLE_PROTEINS
    ]
    selected_example_label = st.selectbox("Примеры белков", example_labels)
    selected_example = EXAMPLE_PROTEINS[example_labels.index(selected_example_label)]

    col_input, col_format = st.columns([1, 1])
    with col_input:
        accession = normalize_accession(
            st.text_input(
                "UniProt ID",
                value=selected_example["id"],
                key=f"accession_{selected_example['id']}",
            )
        )
    with col_format:
        output_label = st.selectbox("Что показать", list(VIEW_OPTIONS.keys()))
        output_format = VIEW_OPTIONS[output_label]

    fields = st.multiselect("Поля для компактной TSV-таблицы", DEFAULT_FIELDS, default=DEFAULT_FIELDS[:7])

    if st.button("Получить данные из UniProt", type="primary"):
        try:
            if output_format == "fasta":
                st.code(fetch_fasta(accession), language="text")
            elif output_format == "tsv":
                st.dataframe(fetch_flat_fields(accession, fields), use_container_width=True)
            else:
                record = fetch_record(accession)
                summary = record_summary(record)
                if output_format == "card":
                    show_summary_card(summary)
                    with st.expander("Последовательность"):
                        st.code(summary["sequence"], language="text")
                elif output_format == "json_tree":
                    st.json(record, expanded=False)
                else:
                    st.code(json.dumps(record, indent=2, ensure_ascii=False), language="json")
        except Exception as error:
            api_error(error)

with tab_search:
    st.header("Поиск в UniProtKB")
    st.write("Эта вкладка нужна, когда точного UniProt ID еще нет и хочется искать по гену, организму, тексту или статусу записи.")

    query = st.text_input("Запрос UniProt", value="gene:INS AND organism_id:9606")
    size = st.slider("Максимум результатов", min_value=1, max_value=50, value=10)
    search_fields = st.multiselect("Какие поля вернуть", DEFAULT_FIELDS, default=DEFAULT_FIELDS[:8], key="search_fields")

    if st.button("Искать в UniProtKB"):
        try:
            text = get_text(
                f"{BASE_URL}/uniprotkb/search",
                params={
                    "query": query,
                    "format": "tsv",
                    "fields": ",".join(search_fields),
                    "size": str(size),
                },
            )
            rows = [line.split("\t") for line in text.strip().splitlines()]
            if len(rows) > 1:
                st.dataframe(pd.DataFrame(rows[1:], columns=rows[0]), use_container_width=True)
            else:
                st.info("Ничего не найдено.")
        except Exception as error:
            api_error(error)

with tab_hits:
    st.header("Карточки похожих белков")
    st.write("Сюда можно вставить ID, которые вернул ваш поиск по эмбеддингам или графу. Формат строки: `UniProtID score`.")

    raw_hits = st.text_area(
        "Кандидаты",
        value="P01308 0.982\nP00698 0.741\nA0A0B4J2F0 0.690",
        height=140,
    )

    if st.button("Собрать карточки"):
        parsed_hits = parse_ids_with_scores(raw_hits)
        if not parsed_hits:
            st.warning("Вставь хотя бы один UniProt ID.")
        for hit in parsed_hits[:20]:
            try:
                record = fetch_record(hit["accession"])
                show_summary_card(record_summary(record), score=hit["score"])
            except Exception as error:
                with st.container(border=True):
                    st.error(f"Не удалось получить {hit['accession']}")
                    st.caption(str(error))

with tab_sequence:
    st.header("Ввод последовательности")
    st.write("Эта вкладка пока не ищет похожие белки. Она готовит вход, который потом можно отправить в ваш сервис эмбеддингов или поиска по графу.")

    sequence_input = st.text_area(
        "Аминокислотная последовательность",
        value="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
        height=180,
    )
    question = st.text_input("Вопрос", value="Что это за белок и какая у него функция?")

    clean_sequence = re.sub(r"[^A-Za-z]", "", sequence_input).upper()
    invalid_chars = sorted(set(re.sub(r"[ACDEFGHIKLMNPQRSTVWYBXZUO]", "", clean_sequence)))

    left, right, third = st.columns(3)
    left.metric("Длина последовательности", len(clean_sequence))
    right.metric("Длина вопроса", len(question))
    third.metric("Неожиданные буквы", len(invalid_chars))

    if invalid_chars:
        st.warning(f"Неожиданные аминокислотные буквы: {', '.join(invalid_chars)}")

    st.code(clean_sequence, language="text")
    st.info(
        "Следующая точка интеграции: отправить очищенную последовательность в сервис эмбеддингов или графового поиска. "
        "Когда он вернет UniProt ID, их можно вставить во вкладку с карточками похожих белков."
    )

with tab_orf:
    st.header("ORF-анализ ДНК/РНК")
    st.write(
        "NCBI ORFfinder доступен как веб-инструмент и standalone-программа, но стабильного публичного REST API у него нет. "
        "Эта вкладка делает локально базовую часть ORFfinder: чистит FASTA, смотрит обе цепи, проверяет 6 рамок считывания и переводит найденные ORF в белок."
    )

    dna_input = st.text_area(
        "Нуклеотидная последовательность",
        value=(
            ">example_insulin_cds\n"
            "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAGGCGCTGTACCTGGTGTGCGGGGAGCGCGGCTTCTTCTACACGCCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGCGCGGGGAGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
        ),
        height=180,
    )
    clean_dna = clean_nucleotide_sequence(dna_input)

    settings_left, settings_mid, settings_right = st.columns(3)
    with settings_left:
        min_nt_length = st.number_input("Минимальная длина ORF, нуклеотиды", min_value=3, max_value=3000, value=90, step=3)
    with settings_mid:
        start_mode = st.selectbox("Где начинать ORF", ["ATG only", "between stops"])
    with settings_right:
        include_partial = st.checkbox("Показывать ORF без стоп-кодона", value=True)
        ignore_nested = st.checkbox("Скрывать вложенные ORF", value=True)

    invalid_nt = sorted(set(re.sub(r"[ACGT]", "", clean_dna)))
    m1, m2, m3 = st.columns(3)
    m1.metric("Длина после очистки", len(clean_dna))
    m2.metric("GC-состав", f"{((clean_dna.count('G') + clean_dna.count('C')) / len(clean_dna) * 100):.1f}%" if clean_dna else "н/д")
    m3.metric("Неожиданные символы", len(invalid_nt))

    if invalid_nt:
        st.warning(f"После очистки остались неожиданные нуклеотидные буквы: {', '.join(invalid_nt)}")

    if st.button("Найти ORF", type="primary"):
        orfs = find_orfs(
            clean_dna,
            min_nt_length=int(min_nt_length),
            start_mode=start_mode,
            include_partial=include_partial,
            ignore_nested=ignore_nested,
        )
        if not orfs:
            st.info("ORF с такими настройками не найдены. Попробуй уменьшить минимальную длину или включить частичные ORF.")
        else:
            table = pd.DataFrame(
                [
                    {
                        "rank": index + 1,
                        "frame": orf["frame"],
                        "strand": orf["strand"],
                        "start": orf["start"],
                        "end": orf["end"],
                        "nt_length": orf["nt_length"],
                        "aa_length": orf["aa_length"],
                        "start_codon": orf["start_codon"],
                        "stop_codon": orf["stop_codon"],
                    }
                    for index, orf in enumerate(orfs[:25])
                ]
            )
            st.dataframe(table, use_container_width=True, hide_index=True)

            best_orf = orfs[0]
            st.subheader("Самый длинный ORF")
            st.caption(
                f"Рамка {best_orf['frame']}, позиции {best_orf['start']}-{best_orf['end']}, "
                f"{best_orf['aa_length']} aa. Этот белок можно отправлять дальше в BLAST, UniProt или embedding-search."
            )
            st.code(f">orf_{best_orf['frame']}_{best_orf['start']}_{best_orf['end']}\n{best_orf['protein']}", language="text")

with tab_docs:
    st.header("Карта API для этого проекта")
    st.markdown(
        f"""
        Основная документация: [{DOCS}]({DOCS})

        Базовый URL API: `{BASE_URL}`

        Полезные endpoint-ы для вашего вероятного сценария:

        | Что нужно | Endpoint | Зачем это нужно |
        |---|---|---|
        | Получить полную запись белка по ID | `/uniprotkb/{{id}}.json` | Превращает найденный embedding-hit в богатую запись о белке. |
        | Получить только последовательность | `/uniprotkb/{{id}}.fasta` | Полезно для alignments, ноутбуков и быстрых проверок. |
        | Искать по гену, организму, ключевому слову, reviewed-статусу | `/uniprotkb/search?query=...` | Нужно, когда точного accession еще нет или надо сравнить кандидатов. |
        | Получить компактную таблицу | `/uniprotkb/search?format=tsv&fields=...` | Быстрее и удобнее для UI, чем полный JSON. |
        | Стримить много результатов | `/uniprotkb/stream?query=...` | Пригодится позже для batch-экспорта или сборки датасета. |
        | Конвертировать ID между базами | `/idmapping/run`, `/idmapping/status/{{jobId}}`, `/idmapping/results/{{jobId}}` | Нужно, если совпадения приходят как RefSeq, Ensembl, PDB, NCBI Gene и т.д. |
        | Смотреть кластеры похожих последовательностей | `/uniref/...` | Полезно, если команде нужен контекст на уровне семейств/кластеров. |
        | Смотреть архив последовательностей | `/uniparc/...` | Полезно для происхождения записей и исторических вариантов. |

        Для первой версии продукта главная роль UniProt такая: ваша система ищет candidate IDs,
        а UniProt объясняет, что это за белки и какие данные вокруг них есть.
        """
    )
