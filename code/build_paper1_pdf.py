#!/usr/bin/env python3
"""
build_paper1_pdf.py — Build paper1_draft_v4.pdf using ReportLab
Fixes §3.2 Table 1 rendering (wide table, 27 rows × 9 columns)

Run: python3 build_paper1_pdf.py
"""

from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.lib.pagesizes import A4, letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageBreak, PageTemplate,
    Paragraph, Spacer, Table, TableStyle, Image,
    KeepTogether, HRFlowable
)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.lib.utils import ImageReader
from reportlab.platypus import Flowable
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


ROOT = Path(__file__).resolve().parent
MD_FILE = ROOT / "paper1_draft_v4.md"
OUT_FILE = ROOT / "paper1_draft_v4.pdf"

# ── Register Unicode-capable TrueType fonts ─────────────────
_ARIAL_PATHS = {
    "Arial":          "/System/Library/Fonts/Supplemental/Arial.ttf",
    "Arial-Bold":     "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "Arial-Italic":   "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
    "Arial-BoldItalic": "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf",
}
# Linux fallback: DejaVu Sans (Unicode-capable, ships with most Linux distros)
_DEJAVU_PATHS = {
    "DejaVu":            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "DejaVu-Bold":       "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "DejaVu-Italic":     "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
    "DejaVu-BoldItalic": "/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf",
}

_USE_TTF = False
_TTF_FAMILY = None  # "Arial" (Mac) or "DejaVu" (Linux) once registered

def _try_register(paths, family_name):
    """Register all four TTFs of a family if all exist; return True on success."""
    if not all(Path(p).exists() for p in paths.values()):
        return False
    try:
        for fname, fpath in paths.items():
            pdfmetrics.registerFont(TTFont(fname, fpath))
        from reportlab.pdfbase.pdfmetrics import registerFontFamily
        registerFontFamily(
            family_name,
            normal=family_name,
            bold=f"{family_name}-Bold",
            italic=f"{family_name}-Italic",
            boldItalic=f"{family_name}-BoldItalic",
        )
        return True
    except Exception:
        return False

if _try_register(_ARIAL_PATHS, "Arial"):
    _USE_TTF = True
    _TTF_FAMILY = "Arial"
elif _try_register(_DEJAVU_PATHS, "DejaVu"):
    _USE_TTF = True
    _TTF_FAMILY = "DejaVu"

if _USE_TTF:
    FONT_NORMAL     = _TTF_FAMILY
    FONT_BOLD       = f"{_TTF_FAMILY}-Bold"
    FONT_ITALIC     = f"{_TTF_FAMILY}-Italic"
    FONT_BOLDITALIC = f"{_TTF_FAMILY}-BoldItalic"
else:
    FONT_NORMAL     = "Helvetica"
    FONT_BOLD       = "Helvetica-Bold"
    FONT_ITALIC     = "Helvetica-Oblique"
    FONT_BOLDITALIC = "Helvetica-BoldOblique"

# ── Page geometry ──────────────────────────────────────────────
PAGE_W, PAGE_H = letter
MARGIN = 2.5 * cm
TEXT_W = PAGE_W - 2 * MARGIN

# ── Styles ──────────────────────────────────────────────────────
styles = getSampleStyleSheet()

TITLE_STYLE = ParagraphStyle(
    "Title",
    parent=styles["Normal"],
    fontSize=14,
    leading=20,
    alignment=TA_CENTER,
    spaceAfter=6,
    fontName=FONT_BOLD,
)
AUTHOR_STYLE = ParagraphStyle(
    "Author",
    parent=styles["Normal"],
    fontSize=11,
    leading=14,
    alignment=TA_CENTER,
    spaceAfter=4,
    fontName=FONT_NORMAL,
)
AFFIL_STYLE = ParagraphStyle(
    "Affil",
    parent=styles["Normal"],
    fontSize=10,
    leading=13,
    alignment=TA_CENTER,
    spaceAfter=12,
    fontName=FONT_ITALIC,
)
H1_STYLE = ParagraphStyle(
    "H1",
    parent=styles["Normal"],
    fontSize=13,
    leading=17,
    spaceBefore=14,
    spaceAfter=4,
    fontName=FONT_BOLD,
    keepWithNext=1,
)
H2_STYLE = ParagraphStyle(
    "H2",
    parent=styles["Normal"],
    fontSize=11,
    leading=15,
    spaceBefore=10,
    spaceAfter=3,
    fontName=FONT_BOLD,
    keepWithNext=1,
)
BODY_STYLE = ParagraphStyle(
    "Body",
    parent=styles["Normal"],
    fontSize=10,
    leading=14,
    alignment=TA_JUSTIFY,
    spaceAfter=6,
    fontName=FONT_NORMAL,
)
ABSTRACT_STYLE = ParagraphStyle(
    "Abstract",
    parent=BODY_STYLE,
    leftIndent=0.5 * cm,
    rightIndent=0.5 * cm,
    fontSize=9.5,
    leading=13,
)
CAPTION_STYLE = ParagraphStyle(
    "Caption",
    parent=BODY_STYLE,
    fontSize=9,
    leading=12,
    fontName=FONT_ITALIC,
    spaceAfter=8,
)
TABLE_HDR_STYLE = ParagraphStyle(
    "TblHdr",
    parent=styles["Normal"],
    fontSize=8,
    leading=10,
    alignment=TA_CENTER,
    fontName=FONT_BOLD,
)
TABLE_CELL_STYLE = ParagraphStyle(
    "TblCell",
    parent=styles["Normal"],
    fontSize=8,
    leading=10,
    alignment=TA_CENTER,
    fontName=FONT_NORMAL,
)
TABLE_CELL_L = ParagraphStyle(
    "TblCellL",
    parent=TABLE_CELL_STYLE,
    alignment=TA_LEFT,
)
TABLE_TITLE_STYLE = ParagraphStyle(
    "TblTitle",
    parent=BODY_STYLE,
    fontName=FONT_BOLD,
    spaceAfter=3,
    spaceBefore=8,
)
REF_STYLE = ParagraphStyle(
    "Ref",
    parent=BODY_STYLE,
    fontSize=9,
    leading=12,
    leftIndent=0.4 * cm,
    firstLineIndent=-0.4 * cm,
    spaceAfter=3,
)


def sanitize_latin(text: str) -> str:
    """Normalize characters for ReportLab. With TTF (Arial): Unicode works,
    just convert special markers to XML markup. Without TTF: ASCII fallbacks."""
    # Always: convert Unicode superscript/subscript markers → ReportLab XML markup
    markup_subs = [
        ('ᵃ', '<super>a</super>'),
        ('ᵇ', '<super>b</super>'),
        ('†', '<super>+</super>'),
        ('₁', '<sub>1</sub>'),
        ('₂', '<sub>2</sub>'),
        ('₀', '<sub>0</sub>'),
        ('¹', '<super>1</super>'),
        ('²', '<super>2</super>'),
        ('³', '<super>3</super>'),
        ('⁴', '<super>4</super>'),
        ('⁵', '<super>5</super>'),
        ('−', '-'),    # proper minus → ASCII hyphen
        ('–', '-'),    # en-dash
        ('—', '--'),   # em-dash
    ]
    for src, dst in markup_subs:
        text = text.replace(src, dst)

    if not _USE_TTF:
        # Helvetica cannot render these; replace with ASCII
        latin_fallbacks = [
            ('ć', 'c'), ('č', 'c'), ('š', 's'), ('ž', 'z'), ('ń', 'n'),
            ('×', 'x'), ('≥', '>='), ('≤', '<='), ('•', '*'),
        ]
        for src, dst in latin_fallbacks:
            text = text.replace(src, dst)
    return text


def md_to_rl(text: str) -> str:
    """Convert inline markdown to ReportLab XML, with character safety."""
    # Bold-italic: ***...***
    text = re.sub(r'\*\*\*(.+?)\*\*\*', r'<b><i>\1</i></b>', text)
    # Bold: **...** or __...__
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'__(.+?)__', r'<b>\1</b>', text)
    # Italic: *...* or _..._
    text = re.sub(r'\*(.+?)\*', r'<i>\1</i>', text)
    text = re.sub(r'_(.+?)_', r'<i>\1</i>', text)
    # Superscript: ^...^
    text = re.sub(r'\^(.+?)\^', r'<super>\1</super>', text)
    # Subscript: ~...~
    text = re.sub(r'~(.+?)~', r'<sub>\1</sub>', text)
    # Escape & for XML (must come before other replacements that insert &)
    text = text.replace('&', '&amp;')
    text = text.replace('&amp;amp;', '&amp;')
    # Sanitize non-Helvetica characters
    text = sanitize_latin(text)
    return text


def clean_cell(text: str) -> str:
    """Strip markdown table cell markup for plain text cells."""
    text = text.strip()
    # Remove surrounding italic markers
    text = re.sub(r'^\*(.+)\*$', r'\1', text)
    # Remove backslash before pipe
    text = text.replace('\\|', '|')
    return text


def build_table1(text_width: float) -> list:
    """
    Build the 27-row Table 1 as a ReportLab flowable list.
    Returns [Paragraph(title), Table, Paragraph(footnotes)]
    """
    title = Paragraph(
        "<b>Table 1.</b> Complete 27-PAH descriptor set with binary "
        "classification at <i>K</i>/|Aut(<i>G</i>)| ≥ 5.0.",
        TABLE_TITLE_STYLE,
    )

    # Column headers (split into two lines where needed)
    headers = [
        Paragraph("#", TABLE_HDR_STYLE),
        Paragraph("Compound", TABLE_HDR_STYLE),
        Paragraph("<i>K</i>", TABLE_HDR_STYLE),
        Paragraph("|Aut(<i>G</i>)|", TABLE_HDR_STYLE),
        Paragraph("<i>K</i>/|Aut|", TABLE_HDR_STYLE),
        Paragraph("PEF", TABLE_HDR_STYLE),
        Paragraph("log<sub>10</sub>\n(PEF)", TABLE_HDR_STYLE),
        Paragraph("Bay/\nFjord", TABLE_HDR_STYLE),
        Paragraph("Classif.<super>a</super>", TABLE_HDR_STYLE),
    ]

    # Data rows: (number, name, K, |Aut|, K/|Aut|, PEF, log10PEF, Bay, Class)
    rows_data = [
        (1,  "Benzene",                      2,  12,  "0.17",  "0.001", "−3.00", "N",  "TN"),
        (2,  "Naphthalene",                  3,   4,  "0.75",  "0.001", "−3.00", "N",  "TN"),
        (3,  "Acenaphthylene",               3,   2,  "1.50",  "0.001", "−3.00", "N",  "TN"),
        (4,  "Fluoranthene",                 6,   2,  "3.00",  "0.001", "−3.00", "N",  "TN"),
        (5,  "Anthracene",                   4,   4,  "1.00",  "0.01",  "−2.00", "N",  "TN"),
        (6,  "Phenanthrene",                 5,   2,  "2.50",  "0.001", "−3.00", "Y",  "TN"),
        (7,  "Pyrene",                       6,   4,  "1.50",  "0.001", "−3.00", "N",  "TN"),
        (8,  "Triphenylene",                 9,   6,  "1.50",  "0.001", "−3.00", "N",  "TN"),
        (9,  "Chrysene",                     8,   2,  "4.00",  "0.01",  "−2.00", "Y",  "TN"),
        (10, "Benz[a]anthracene",            7,   1,  "7.00",  "0.1",   "−1.00", "Y",  "TP"),
        (11, "Benzo[c]phenanthrene",         8,   2,  "4.00",  "0.023", "−1.64", "Y†", "TN"),
        (12, "Benzo[a]pyrene",               9,   1,  "9.00",  "1.0",   " 0.00", "Y",  "TP"),
        (13, "3-Methylcholanthrene",         7,   1,  "7.00",  "1.833", "+0.26", "Y",  "TP"),
        (14, "5-Methylchrysene",             8,   2,  "4.00",  "1.0",   " 0.00", "Y",  "FN"),
        (15, "DMBA",                         7,   1,  "7.00",  "20.83", "+1.32", "Y",  "TP"),
        (16, "Dibenz[a,h]anthracene",       10,   2,  "5.00",  "5.0",   "+0.70", "Y",  "TP"),
        (17, "Dibenzo[a,l]pyrene",          16,   1, "16.00",  "10.0",  "+1.00", "Y",  "TP"),
        (18, "Benzo[b]fluoranthene",        10,   1, "10.00",  "0.1",   "−1.00", "Y",  "TP"),
        (19, "Benzo[k]fluoranthene",         9,   2,  "4.50",  "0.1",   "−1.00", "Y",  "FN"),
        (20, "Perylene",                     9,   4,  "2.25",  "0.001", "−3.00", "N",  "TN"),
        (21, "Benzo[ghi]perylene",          14,   2,  "7.00",  "0.01",  "−2.00", "N",  "FP"),
        (22, "Coronene",                    20,  12,  "1.67",  "0.001", "−3.00", "N",  "TN"),
        (23, "Indeno[1,2,3-cd]pyrene",      12,   1, "12.00",  "0.1",   "−1.00", "Y",  "TP"),
        (24, "Dibenzo[a,e]pyrene",          17,   1, "17.00",  "1.0",   " 0.00", "Y",  "TP"),
        (25, "Dibenzo[a,h]pyrene",          13,   2,  "6.50",  "10.0",  "+1.00", "Y",  "TP"),
        (26, "Dibenzo[a,i]pyrene",          14,   2,  "7.00",  "10.0",  "+1.00", "Y",  "TP"),
        (27, "Benzo[e]pyrene",              11,   2,  "5.50",  "0.001", "−3.00", "N",  "FP"),
    ]

    # Classification color coding
    CLASSIF_COLORS = {
        "TP": colors.HexColor("#d4edda"),   # light green
        "TN": colors.HexColor("#d1ecf1"),   # light blue
        "FP": colors.HexColor("#ffeeba"),   # light yellow
        "FN": colors.HexColor("#f8d7da"),   # light pink/red
    }

    def make_cell(val, style=TABLE_CELL_STYLE):
        return Paragraph(str(val), style)

    # Build table data
    table_data = [headers]
    classif_col = []
    for row in rows_data:
        num, name, K, aut, k_aut, pef, log_pef, bay, classif = row
        classif_col.append(classif)
        table_data.append([
            make_cell(num),
            make_cell(name, TABLE_CELL_L),
            make_cell(K),
            make_cell(aut),
            make_cell(k_aut),
            make_cell(pef),
            make_cell(log_pef),
            make_cell(bay),
            make_cell(classif),
        ])

    # Column widths (total must = text_width)
    # #(0.3) | Compound(3.8cm) | K(0.6) | Aut(0.8) | K/Aut(0.9) | PEF(0.9) | logPEF(0.9) | Bay(0.7) | Class(0.9)
    # Let's size in cm, then scale to fit text_width
    col_w_raw = [0.30, 3.80, 0.55, 0.70, 0.80, 0.80, 0.80, 0.60, 0.75]
    total_raw = sum(col_w_raw)  # in cm
    col_widths = [w * cm * (text_width / (total_raw * cm)) for w in col_w_raw]

    tbl = Table(table_data, colWidths=col_widths, repeatRows=1)

    # Base style
    ts = [
        ("BACKGROUND",   (0, 0), (-1, 0),  colors.HexColor("#343a40")),
        ("TEXTCOLOR",    (0, 0), (-1, 0),  colors.white),
        ("ALIGN",        (0, 0), (-1, -1), "CENTER"),
        ("ALIGN",        (1, 0), (1, -1),  "LEFT"),
        ("FONTNAME",     (0, 0), (-1, 0),  FONT_BOLD),
        ("FONTSIZE",     (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8f9fa")]),
        ("GRID",         (0, 0), (-1, -1), 0.3, colors.HexColor("#dee2e6")),
        ("TOPPADDING",   (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING",(0, 0), (-1, -1), 3),
        ("LEFTPADDING",  (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("VALIGN",       (0, 0), (-1, -1), "MIDDLE"),
    ]

    # Color-code classification column (index 8)
    for i, classif in enumerate(classif_col):
        row_idx = i + 1  # 0 = header
        bg = CLASSIF_COLORS.get(classif[:2], colors.white)
        ts.append(("BACKGROUND", (8, row_idx), (8, row_idx), bg))

    tbl.setStyle(TableStyle(ts))

    # Footnotes are in the markdown source and will be picked up by the
    # standalone footnote handler in parse_markdown(). Do NOT add them here
    # to avoid duplication.
    return [title, tbl]


def add_figure(path_str: str, caption_text: str, max_width: float, max_height: float) -> list:
    """Return [Image, Paragraph(caption)] or [] if file not found."""
    # Strip URL encoding and markdown image syntax
    path_str = path_str.replace("%20", " ")
    p = Path(path_str)
    if not p.exists():
        return [Paragraph(f"[Figure not found: {p.name}]", CAPTION_STYLE)]
    try:
        img = ImageReader(str(p))
        iw, ih = img.getSize()
        if iw and ih:
            scale = min(max_width / iw, max_height / ih, 1.0)
            w, h = iw * scale, ih * scale
        else:
            w, h = max_width, max_height * 0.6
        flowable_img = Image(str(p), width=w, height=h)
        flowable_img.hAlign = "CENTER"
    except Exception as e:
        return [Paragraph(f"[Image error: {e}]", CAPTION_STYLE)]

    caption = Paragraph(md_to_rl(caption_text), CAPTION_STYLE)
    return [Spacer(1, 4), flowable_img, Spacer(1, 4), caption]


def parse_markdown(md_text: str, text_width: float) -> list:
    """Parse markdown and produce a list of ReportLab flowables."""
    story = []
    lines = md_text.split("\n")
    i = 0
    in_abstract = False
    abstract_lines = []
    in_references = False
    ref_counter = 0

    def flush_abstract():
        nonlocal abstract_lines, in_abstract
        if abstract_lines:
            text = " ".join(abstract_lines).strip()
            story.append(Paragraph(md_to_rl(text), ABSTRACT_STYLE))
            story.append(Spacer(1, 6))
            abstract_lines = []
        in_abstract = False

    while i < len(lines):
        line = lines[i]

        # ── Title (first H1 in file) ──────────────────────────
        if line.startswith("# ") and not story:
            title_text = line[2:].strip()
            story.append(Spacer(1, 0.3 * cm))
            story.append(Paragraph(md_to_rl(title_text), TITLE_STYLE))
            i += 1
            continue

        # ── Author / affiliation ──────────────────────────────
        if i == 2 and line.strip() and not line.startswith("#"):
            story.append(Paragraph(line.strip(), AUTHOR_STYLE))
            i += 1
            continue
        if i == 4 and line.strip() and not line.startswith("#"):
            story.append(Paragraph(line.strip(), AFFIL_STYLE))
            story.append(HRFlowable(width="100%", thickness=0.5, color=colors.grey))
            story.append(Spacer(1, 6))
            i += 1
            continue

        # ── HR ────────────────────────────────────────────────
        if line.strip() == "---":
            if in_abstract:
                flush_abstract()
            story.append(HRFlowable(width="100%", thickness=0.5, color=colors.grey))
            story.append(Spacer(1, 4))
            i += 1
            continue

        # ── Section headings ─────────────────────────────────
        if line.startswith("## "):
            if in_abstract:
                flush_abstract()
            heading = line[3:].strip()
            if heading == "Abstract":
                story.append(Paragraph("Abstract", H1_STYLE))
                in_abstract = True
                i += 1
                continue
            if heading == "References":
                in_references = True
            story.append(Paragraph(md_to_rl(heading), H1_STYLE))
            i += 1
            continue

        if line.startswith("### "):
            if in_abstract:
                flush_abstract()
            story.append(Paragraph(md_to_rl(line[4:].strip()), H2_STYLE))
            i += 1
            continue

        # ── Figure image ─────────────────────────────────────
        if line.strip().startswith("!["):
            if in_abstract:
                flush_abstract()
            # Extract alt text and path
            m = re.match(r'!\[([^\]]*)\]\(([^)]+)\)', line.strip())
            if m:
                # Next non-empty line is the caption
                caption = ""
                j = i + 1
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines) and not lines[j].startswith("!") and not lines[j].startswith("#") and not lines[j].startswith("|"):
                    caption = lines[j].strip()
                    i = j + 1
                else:
                    caption = m.group(1)
                    i += 1
                story.extend(add_figure(m.group(2), caption, text_width, 4.5 * inch))
            else:
                i += 1
            continue

        # ── Placeholder figures [Figure N: ...] ─────────────
        if re.match(r'\*\*\[Figure \d', line.strip()):
            i += 1
            continue

        # ── Table (pipe-delimited markdown) ──────────────────
        if line.strip().startswith("|"):
            if in_abstract:
                flush_abstract()
            # Collect all lines of this table
            tbl_lines = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                tbl_lines.append(lines[i])
                i += 1
            # Check if this follows a "Table 1." bold line (already handled)
            # Check for separator row (---|---)
            sep_idx = None
            for ti, tl in enumerate(tbl_lines):
                if re.match(r'[\|\s\-:]+$', tl.replace('|', '').replace('-', '').replace(':', '').replace(' ', '') or 'x'):
                    pass
                if all(c in '-| :' for c in tl):
                    sep_idx = ti
                    break

            if sep_idx is None:
                i_start, i_data = 0, 0
            else:
                i_start = 0
                i_data = sep_idx + 1

            def parse_row(row_str):
                cells = [c.strip() for c in row_str.split("|")]
                # Remove empty first/last if pipe-delimited
                if cells and not cells[0]:
                    cells = cells[1:]
                if cells and not cells[-1]:
                    cells = cells[:-1]
                return cells

            header_cells = parse_row(tbl_lines[i_start])
            data_rows = [parse_row(tbl_lines[j]) for j in range(i_data, len(tbl_lines))]

            n_col = len(header_cells)
            if n_col == 0:
                continue

            # Heuristic: the 27-PAH compound table (Table 1) is the ONLY
            # wide table pre-rendered via build_table1. It is identified by
            # a first-column header of "#". All other pipe-tables (including
            # the §3.2 Bonferroni×17 table with 7 columns) should render inline.
            first_hdr = (header_cells[0] if header_cells else "").strip()
            is_compound_table = first_hdr in {"#", "No.", "No"}

            if is_compound_table:
                # Skip — Table 1 already appended via build_table1 trigger
                pass
            else:
                # Build small table
                col_w = text_width / n_col

                def rl_cells(cells, is_hdr=False):
                    sty = TABLE_HDR_STYLE if is_hdr else TABLE_CELL_STYLE
                    return [Paragraph(md_to_rl(clean_cell(c)), sty) for c in cells]

                t_data = [rl_cells(header_cells, True)]
                for dr in data_rows:
                    # Pad or trim to n_col
                    dr = (dr + [""] * n_col)[:n_col]
                    t_data.append(rl_cells(dr))

                small_tbl = Table(t_data, colWidths=[col_w] * n_col, repeatRows=1)
                small_tbl.setStyle(TableStyle([
                    ("BACKGROUND",   (0, 0), (-1, 0), colors.HexColor("#343a40")),
                    ("TEXTCOLOR",    (0, 0), (-1, 0), colors.white),
                    ("GRID",         (0, 0), (-1, -1), 0.4, colors.HexColor("#dee2e6")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                     [colors.white, colors.HexColor("#f8f9fa")]),
                    ("TOPPADDING",   (0, 0), (-1, -1), 3),
                    ("BOTTOMPADDING",(0, 0), (-1, -1), 3),
                    ("LEFTPADDING",  (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("FONTSIZE",     (0, 0), (-1, -1), 9),
                    ("VALIGN",       (0, 0), (-1, -1), "MIDDLE"),
                    ("ALIGN",        (0, 0), (-1, -1), "CENTER"),
                ]))
                story.append(small_tbl)
                story.append(Spacer(1, 4))
            continue

        # ── Footnote lines (start with ᵃ † ¹ ²) ─────────────
        if line.strip() and line.strip()[0] in ('ᵃ', '†', '¹', '²', '³', '⁴', '⁵'):
            if in_abstract:
                flush_abstract()
            story.append(Paragraph(md_to_rl(line.strip()), CAPTION_STYLE))
            i += 1
            continue

        # ── Bold "Table N." line ──────────────────────────────
        if re.match(r'\*\*Table \d', line.strip()):
            if in_abstract:
                flush_abstract()
            # Insert hard-coded Table 1 here
            story.extend(build_table1(text_width))
            story.append(Spacer(1, 6))
            i += 1
            continue

        # ── Ordered list items ────────────────────────────────
        if re.match(r'^\d+\. ', line.strip()):
            if in_abstract:
                flush_abstract()
            if in_references:
                story.append(Paragraph(md_to_rl(line.strip()), REF_STYLE))
            else:
                story.append(Paragraph(md_to_rl(line.strip()), BODY_STYLE))
            i += 1
            continue

        # ── Unordered list items ──────────────────────────────
        if line.strip().startswith("- ") or line.strip().startswith("* "):
            if in_abstract:
                flush_abstract()
            item = line.strip()[2:]
            story.append(Paragraph("• " + md_to_rl(item),
                                   ParagraphStyle("li", parent=BODY_STYLE,
                                                  leftIndent=0.4*cm)))
            i += 1
            continue

        # ── Empty line ────────────────────────────────────────
        if not line.strip():
            if in_abstract and abstract_lines:
                # paragraph break within abstract
                abstract_lines.append(" ")
            i += 1
            continue

        # ── Regular paragraph ─────────────────────────────────
        if in_abstract:
            abstract_lines.append(line.strip())
        else:
            story.append(Paragraph(md_to_rl(line.strip()), BODY_STYLE))

        i += 1

    if in_abstract:
        flush_abstract()

    return story


def build_pdf():
    md_text = MD_FILE.read_text(encoding="utf-8")

    story = parse_markdown(md_text, TEXT_W)

    # ── Document ──────────────────────────────────────────────
    doc = BaseDocTemplate(
        str(OUT_FILE),
        pagesize=letter,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN,
    )

    frame = Frame(
        MARGIN, MARGIN,
        PAGE_W - 2 * MARGIN, PAGE_H - 2 * MARGIN,
        id="main",
    )

    def on_page(canvas, doc):
        canvas.saveState()
        canvas.setFont(FONT_NORMAL, 8)
        canvas.setFillColor(colors.grey)
        canvas.drawCentredString(PAGE_W / 2, 0.8 * cm, str(doc.page))
        canvas.restoreState()

    template = PageTemplate(id="main", frames=[frame], onPage=on_page)
    doc.addPageTemplates([template])

    doc.build(story)
    size_kb = OUT_FILE.stat().st_size // 1024
    print(f"✓ Built: {OUT_FILE}")
    print(f"  Size: {size_kb} KB")
    print(f"  Pages: approximately {len(story) // 20 + 1}")


if __name__ == "__main__":
    build_pdf()
