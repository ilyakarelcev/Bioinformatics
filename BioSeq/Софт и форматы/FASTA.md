# FASTA

**FASTA** — это простой текстовый формат для записи биологических последовательностей: ДНК, РНК или белков. Первая строка начинается с `>` и содержит описание, а следующие строки содержат саму последовательность.

Пример: инсулин человека из [[_UniProt|UniProt]].

```fasta
>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT
RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

## Разбор заголовка

| Часть заголовка | Что означает |
|---|---|
| `>` | Начало FASTA-заголовка. |
| `sp` | UniProtKB/[[Swiss-Prot]]: вручную проверенная запись. |
| `P01308` | [[Accession number]]: стабильный идентификатор записи в UniProt. |
| `INS_HUMAN` | Entry name: короткое имя записи; `INS` — инсулин, `HUMAN` — человек. |
| `Insulin` | Название белка. |
| `OS=Homo sapiens` | Organism Species: организм, из которого получена последовательность. |
| `OX=9606` | Organism Taxonomy ID: таксономический номер человека. |
| `GN=INS` | Gene Name: название гена. |
| `PE=1` | [[Protein Existence]]: доказательность существования белка; `1` — подтвержден на уровне белка. |
| `SV=1` | Sequence Version: версия последовательности. |

Вариант заголовка из [[NCBI]] может быть короче:

```fasta
>NP_000198.1 insulin [Homo sapiens]
MALWMRLLPLLALLALWGPDPAAA
```
