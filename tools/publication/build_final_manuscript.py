#!/usr/bin/env python3
"""Build the G-Synth manuscript and Supporting Information DOCX files.

The source insulin manuscript is evidence, not a layout template.  Documents use
the standard_business_brief preset with a named "journal manuscript" override:
single-column US Letter, restrained colour, scientific captions and compact
references.  Values are applied explicitly instead of relying on Word defaults.
"""

from __future__ import annotations

import json
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[2]
PUBLICATION = ROOT / "publication"
EVIDENCE = ROOT / "publication_evidence"
FIG = EVIDENCE / "manuscript_figures"
DESIGN_JSON = EVIDENCE / "glargine_ab_design_and_cloning.json"
SEQ_JSON = EVIDENCE / "sequencing_validation" / "glargine_approved_trace_validation.json"
PRIMARY_RECORD_JSON = ROOT / "tools" / "publication" / "primary_acs_draft_record.json"
SOURCE_FIG = EVIDENCE / "source_experimental_figures"
INTERFACE = EVIDENCE / "interface_evidence"

MAIN_OUT = PUBLICATION / "G-Synth_Final_Master_Manuscript.docx"
SI_OUT = PUBLICATION / "G-Synth_Master_Supporting_Information.docx"

NAVY = "12233F"
TEAL = "39757D"
PURPLE = "705493"
GREEN = "4C8B66"
MUTED = "5D6B78"
LIGHT = "F2F4F7"
PALE_TEAL = "E8F1F3"
PALE_GREEN = "EAF3EC"
CAUTION = "FFF4CE"
WHITE = "FFFFFF"
BLACK = "111111"


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=120, bottom=80, end=120) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for tag, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{tag}"))
        if node is None:
            node = OxmlElement(f"w:{tag}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_table_borders(table, color="C8D0D8", size="6") -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.find(qn("w:tblBorders"))
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        el = borders.find(qn(f"w:{edge}"))
        if el is None:
            el = OxmlElement(f"w:{edge}")
            borders.append(el)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), size)
        el.set(qn("w:color"), color)


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    tr_pr.append(repeat)


def cant_split(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:cantSplit")
    tr_pr.append(el)


def set_width(cell, width_in: float) -> None:
    cell.width = Inches(width_in)
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(round(width_in * 1440)))
    tc_w.set(qn("w:type"), "dxa")


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("")
    run.font.name = "Calibri"
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor.from_string(MUTED)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    paragraph._p.append(fld)


def set_keep_with_next(paragraph, value=True) -> None:
    paragraph.paragraph_format.keep_with_next = value


def configure_document(doc: Document, running_title: str) -> None:
    properties = doc.core_properties
    properties.title = running_title
    properties.subject = "G-Synth scientific software publication"
    properties.author = "Mohamed Merzoug"
    properties.last_modified_by = "Mohamed Merzoug"
    properties.keywords = "G-Synth, synthetic biology, nucleic-acid design, software validation"
    properties.revision = 1

    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.85)
    section.bottom_margin = Inches(0.85)
    section.left_margin = Inches(0.9)
    section.right_margin = Inches(0.9)
    section.header_distance = Inches(0.42)
    section.footer_distance = Inches(0.42)

    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Calibri"
    normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor.from_string(BLACK)
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.10

    for name, size, color, before, after in (
        ("Heading 1", 16, NAVY, 16, 8),
        ("Heading 2", 13, TEAL, 12, 6),
        ("Heading 3", 11.5, PURPLE, 8, 4),
    ):
        s = styles[name]
        s.font.name = "Calibri"
        s._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
        s._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
        s.font.size = Pt(size)
        s.font.bold = True
        s.font.color.rgb = RGBColor.from_string(color)
        s.paragraph_format.space_before = Pt(before)
        s.paragraph_format.space_after = Pt(after)
        s.paragraph_format.keep_with_next = True

    cap = styles["Caption"]
    cap.font.name = "Calibri"
    cap.font.size = Pt(9)
    cap.font.color.rgb = RGBColor.from_string(BLACK)
    cap.font.italic = False
    cap.paragraph_format.space_before = Pt(3)
    cap.paragraph_format.space_after = Pt(8)
    cap.paragraph_format.line_spacing = 1.0
    cap.paragraph_format.keep_with_next = False

    header = section.header.paragraphs[0]
    header.text = running_title
    header.alignment = WD_ALIGN_PARAGRAPH.LEFT
    header.paragraph_format.space_after = Pt(2)
    for run in header.runs:
        run.font.name = "Calibri"
        run.font.size = Pt(8.5)
        run.font.bold = True
        run.font.color.rgb = RGBColor.from_string(MUTED)
    p_pr = header._p.get_or_add_pPr()
    p_bdr = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), "6")
    bottom.set(qn("w:space"), "3")
    bottom.set(qn("w:color"), "D9E0E6")
    p_bdr.append(bottom)
    p_pr.append(p_bdr)
    add_page_number(section.footer.paragraphs[0])


def add_title_block(doc: Document, title: str, subtitle: str | None = None) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_before = Pt(14)
    p.paragraph_format.space_after = Pt(5)
    r = p.add_run(title)
    r.bold = True
    r.font.name = "Calibri"
    r.font.size = Pt(20)
    r.font.color.rgb = RGBColor.from_string(NAVY)
    if subtitle:
        p2 = doc.add_paragraph()
        p2.paragraph_format.space_after = Pt(12)
        r2 = p2.add_run(subtitle)
        r2.font.name = "Calibri"
        r2.font.size = Pt(12)
        r2.font.color.rgb = RGBColor.from_string(TEAL)
        r2.italic = True


def add_authors(doc: Document) -> None:
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(5)
    r = p.add_run(
        "Mohamed Merzougᵃ,* · Zohra Yasmine Zaterᵃ,ᵇ · Yasmine Saidiᵃ · "
        "Amaria Ilhem Hammadiᵃ · Marwa Airecheᵃ · Keltoum Bendidaᵃ · "
        "Hadjer Soumia Bouderbalaᵃ · Soheir Bouzidiᵃ · Djamal Saidiᵃ"
    )
    r.bold = True
    r.font.size = Pt(10.5)
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    p.add_run(
        "ᵃ Higher School of Biological Sciences of Oran, BP 1042 Saim Mohamed, "
        "Cité Emir Abdelkader, 31000 Oran, Algeria"
    ).font.size = Pt(9.5)
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    p.add_run(
        "ᵇ Oran 1 University Ahmed Ben Bella, BP 1524 El M'Naouer, 31000 Oran, Algeria"
    ).font.size = Pt(9.5)
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(12)
    rr = p.add_run("* Corresponding author: Mohamed Merzoug; merzoug.mohamed@essb-oran.edu.dz")
    rr.font.size = Pt(9.5)
    rr.italic = True


def add_text(doc: Document, text: str, bold_lead: str | None = None, italic=False) -> None:
    p = doc.add_paragraph()
    if bold_lead and text.startswith(bold_lead):
        p.add_run(bold_lead).bold = True
        p.add_run(text[len(bold_lead):])
    else:
        p.add_run(text)
    if italic:
        for r in p.runs:
            r.italic = True


def add_bullets(doc: Document, items: list[str]) -> None:
    for item in items:
        p = doc.add_paragraph(style="List Bullet")
        p.paragraph_format.left_indent = Inches(0.5)
        p.paragraph_format.first_line_indent = Inches(-0.25)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing = 1.10
        p.add_run(item)


def add_callout(doc: Document, title: str, body: str, fill=PALE_TEAL) -> None:
    table = doc.add_table(rows=1, cols=1)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_repeat_table_header(table.rows[0])
    cant_split(table.rows[0])
    cell = table.cell(0, 0)
    set_width(cell, 6.6)
    set_cell_margins(cell, top=120, bottom=120, start=150, end=150)
    set_cell_shading(cell, fill)
    set_table_borders(table, color="B8C8CF", size="6")
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(3)
    r = p.add_run(title)
    r.bold = True
    r.font.color.rgb = RGBColor.from_string(NAVY)
    p2 = cell.add_paragraph(body)
    p2.paragraph_format.space_after = Pt(0)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)


def add_table(doc: Document, headers: list[str], rows: list[list[str]], widths: list[float], font_size=8.7) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    set_table_borders(table)
    hdr = table.rows[0]
    set_repeat_table_header(hdr)
    for cell, label, width in zip(hdr.cells, headers, widths, strict=True):
        set_width(cell, width)
        set_cell_margins(cell)
        set_cell_shading(cell, NAVY)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(label)
        r.bold = True
        r.font.size = Pt(font_size)
        r.font.color.rgb = RGBColor.from_string(WHITE)
    for ri, data in enumerate(rows):
        row = table.add_row()
        cant_split(row)
        for cell, value, width in zip(row.cells, data, widths, strict=True):
            set_width(cell, width)
            set_cell_margins(cell)
            if ri % 2:
                set_cell_shading(cell, LIGHT)
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(str(value))
            r.font.size = Pt(font_size)
    doc.add_paragraph().paragraph_format.space_after = Pt(1)


def add_figure(doc: Document, path: Path, caption: str, width=6.6) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_with_next = True
    picture = p.add_run().add_picture(str(path), width=Inches(width))
    picture._inline.docPr.set("title", path.stem)
    picture._inline.docPr.set("descr", caption)
    cap = doc.add_paragraph(style="Caption")
    cap.add_run(caption)


def add_reference_list(doc: Document, references: list[str]) -> None:
    for i, ref in enumerate(references, 1):
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.25)
        p.paragraph_format.first_line_indent = Inches(-0.25)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(f"{i}. {ref}")
        r.font.size = Pt(8.5)


def main_references(primary: dict) -> list[str]:
    """Retain every reference from the primary draft, then add software literature."""
    return primary["references"] + [
        "Hillson, N. J.; Rosengarten, R. D.; Keasling, J. D. j5 DNA Assembly Design Automation Software. ACS Synth. Biol. 2012, 1, 14–21. https://doi.org/10.1021/sb2000116.",
        "Haines, M. C.; Carling, B.; Marshall, J.; et al. basicsynbio and the BASIC SEVA Collection: Software and Vectors for an Established DNA Assembly Method. Synth. Biol. 2022, 7, ysac023. https://doi.org/10.1093/synbio/ysac023.",
        "Storch, M.; Haines, M. C.; Baldwin, G. S. DNA-BOT: A Low-Cost, Automated DNA Assembly Platform for Synthetic Biology. Synth. Biol. 2020, 5, ysaa010. https://doi.org/10.1093/synbio/ysaa010.",
        "Richardson, S. M.; Wheelan, S. J.; Yarrington, R. M.; Boeke, J. D. GeneDesign: Rapid, Automated Design of Multikilobase Synthetic Genes. Genome Res. 2006, 16, 550–556. https://doi.org/10.1101/gr.4431306.",
        "Richardson, S. M.; Nunley, P. W.; Yarrington, R. M.; Boeke, J. D.; Bader, J. S. GeneDesign 3.0 Is an Updated Synthetic Biology Toolkit. Nucleic Acids Res. 2010, 38, 2603–2606. https://doi.org/10.1093/nar/gkq143.",
        "Zulkower, V.; Rosser, S. DNA Chisel, a Versatile Sequence Optimizer. Bioinformatics 2020, 36, 4508–4509. https://doi.org/10.1093/bioinformatics/btaa558.",
        "Ham, T. S.; Dmytriv, Z.; Plahar, H.; Chen, J.; Hillson, N. J.; Keasling, J. D. Design, Implementation and Practice of JBEI-ICE: An Open Source Biological Part Registry Platform and Tools. Nucleic Acids Res. 2012, 40, e141. https://doi.org/10.1093/nar/gks531.",
        "Kearse, M.; Moir, R.; Wilson, A.; et al. Geneious Basic: An Integrated and Extendable Desktop Software Platform for the Organization and Analysis of Sequence Data. Bioinformatics 2012, 28, 1647–1649. https://doi.org/10.1093/bioinformatics/bts199.",
        "Cock, P. J. A.; Antao, T.; Chang, J. T.; et al. Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics. Bioinformatics 2009, 25, 1422–1423. https://doi.org/10.1093/bioinformatics/btp163.",
        "Galdzicki, M.; Clancy, K. P.; Oberortner, E.; et al. The Synthetic Biology Open Language Provides a Community Standard for Communicating Designs in Synthetic Biology. Nat. Biotechnol. 2014, 32, 545–550. https://doi.org/10.1038/nbt.2891.",
        "Madsen, C.; Goñi-Moreno, A.; Palchick, Z.; et al. Synthetic Biology Open Language Visual Version 2.2. J. Integr. Bioinform. 2019, 16, 20180101. https://doi.org/10.1515/jib-2018-0101.",
        "Plahar, H. A.; et al. Vision and Development of a Design, Implementation, and Verification Automation (DIVA) Software Platform for DNA Construction. ACS Synth. Biol. 2026, 15, 3498–3503. https://doi.org/10.1021/acssynbio.6c00197.",
        "Jeon, J.; et al. CloneCoordinate: A Computational Platform for Coordinating DNA Construction. ACS Synth. Biol. 2025, 14, 4802–4818. https://doi.org/10.1021/acssynbio.5c00582.",
        "Galez, H.; et al. InSillyClo: Software-Assisted Planning of Golden Gate and MoClo Workflows. ACS Synth. Biol. 2026, 15, 353–358. https://doi.org/10.1021/acssynbio.5c00553.",
        "Barker, M.; Chue Hong, N. P.; Katz, D. S.; et al. Introducing the FAIR Principles for Research Software. Sci. Data 2022, 9, 622. https://doi.org/10.1038/s41597-022-01710-x.",
        "Taschuk, M.; Wilson, G. Ten Simple Rules for Making Research Software More Robust. PLoS Comput. Biol. 2017, 13, e1005412. https://doi.org/10.1371/journal.pcbi.1005412.",
        "Brack, P.; Crowther, P.; Soiland-Reyes, S.; et al. Ten Simple Rules for Making a Software Tool Workflow-Ready. PLoS Comput. Biol. 2022, 18, e1009823. https://doi.org/10.1371/journal.pcbi.1009823.",
        "Jiménez, R. C.; Kuzak, M.; Alhamdoosh, M.; et al. Four Simple Recommendations to Encourage Best Practices in Research Software. F1000Research 2017, 6, 876. https://doi.org/10.12688/f1000research.11407.1.",
        "Smith, A. M.; Katz, D. S.; Niemeyer, K. E.; FORCE11 Software Citation Working Group. Software Citation Principles. PeerJ Comput. Sci. 2016, 2, e86. https://doi.org/10.7717/peerj-cs.86.",
        "Sanger, F.; Nicklen, S.; Coulson, A. R. DNA Sequencing with Chain-Terminating Inhibitors. Proc. Natl. Acad. Sci. U.S.A. 1977, 74, 5463–5467. https://doi.org/10.1073/pnas.74.12.5463.",
        "Ewing, B.; Hillier, L.; Wendl, M. C.; Green, P. Base-Calling of Automated Sequencer Traces Using Phred. I. Accuracy Assessment. Genome Res. 1998, 8, 175–185. https://doi.org/10.1101/gr.8.3.175.",
        "Ewing, B.; Green, P. Base-Calling of Automated Sequencer Traces Using Phred. II. Error Probabilities. Genome Res. 1998, 8, 186–194. https://doi.org/10.1101/gr.8.3.186.",
        "SantaLucia, J., Jr. A Unified View of Polymer, Dumbbell, and Oligonucleotide DNA Nearest-Neighbor Thermodynamics. Proc. Natl. Acad. Sci. U.S.A. 1998, 95, 1460–1465. https://doi.org/10.1073/pnas.95.4.1460.",
        "Nakamura, Y.; Gojobori, T.; Ikemura, T. Codon Usage Tabulated from the International DNA Sequence Databases: Status for the Year 2000. Nucleic Acids Res. 2000, 28, 292. https://doi.org/10.1093/nar/28.1.292.",
        "Sharp, P. M.; Li, W.-H. The Codon Adaptation Index—A Measure of Directional Synonymous Codon Usage Bias, and Its Potential Applications. Nucleic Acids Res. 1987, 15, 1281–1295. https://doi.org/10.1093/nar/15.3.1281.",
        "Alexaki, A.; Kames, J.; Holcomb, D. D.; et al. Codon and Codon-Pair Usage Tables (CoCoPUTs): Facilitating Genetic Variation Analyses and Recombinant Gene Design. J. Mol. Biol. 2019, 431, 2434–2441. https://doi.org/10.1016/j.jmb.2019.04.021.",
        "Athey, J.; Alexaki, A.; Osipova, E.; et al. A New and Updated Resource for Codon Usage Tables. BMC Bioinformatics 2017, 18, 391. https://doi.org/10.1186/s12859-017-1793-7.",
        "Ranaghan, M. J.; Li, J. J.; Laprise, D. M.; et al. Assessing Optimal: Inequalities in Codon Optimization Algorithms. BMC Biol. 2021, 19, 36. https://doi.org/10.1186/s12915-021-00968-8.",
        "Subramanian, K.; Payne, B.; Feyertag, F.; Alvarez-Ponce, D. The Codon Statistics Database: A Database of Codon Usage Bias. Mol. Biol. Evol. 2022, 39, msac157. https://doi.org/10.1093/molbev/msac157.",
        "O'Leary, N. A.; Wright, M. W.; Brister, J. R.; et al. Reference Sequence (RefSeq) Database at NCBI: Current Status, Taxonomic Expansion, and Functional Annotation. Nucleic Acids Res. 2016, 44, D733–D745. https://doi.org/10.1093/nar/gkv1189.",
        "National Center for Biotechnology Information. The Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-09-02).",
        "Ben-Bassat, A.; Bauer, K.; Chang, S. Y.; et al. Processing of the Initiation Methionine from Proteins: Properties of the Escherichia coli Methionine Aminopeptidase and Its Gene Structure. J. Bacteriol. 1987, 169, 751–757. https://doi.org/10.1128/jb.169.2.751-757.1987.",
        "Xiao, Q.; Zhang, F.; Nacev, B. A.; Liu, J. O.; Pei, D. Protein N-Terminal Processing: Substrate Specificity of Escherichia coli and Human Methionine Aminopeptidases. Biochemistry 2010, 49, 5588–5599. https://doi.org/10.1021/bi1005464.",
    ]


def build_main(design: dict, seq: dict, primary: dict) -> None:
    doc = Document()
    biopython_version = seq.get("analysis", {}).get("biopython_version", "version recorded in the evidence JSON")
    approved = seq["chains"]
    a_consensus = approved["A"]["gsynth_forward_reverse_consensus_raw_q0"]
    b_consensus = approved["B"]["gsynth_forward_reverse_consensus_raw_q0"]
    configure_document(doc, "G-Synth — end-to-end nucleic-acid design and validation")
    add_title_block(
        doc,
        "G-Synth: Small Sequence Design and Extended Sequence Design for auditable synthesis-ready DNA and post-sequencing validation",
        "Software implementation and experimental validation with PCR-free insulin glargine A/B-chain constructs",
    )
    add_authors(doc)
    doc.add_heading("Abstract", level=1)
    add_text(doc,
        f"Nucleic-acid design for synthesis is often fragmented across codon tools, restriction maps, plasmid editors and sequencing viewers, complicating proof that an ordered molecule is the one later cloned and assessed. We developed G-Synth, an open-source application and Python engine implementing 15 source-traceable host codon profiles and two synthesis-centred workflows: Small Sequence Design (SSD) for one complementary oligonucleotide pair and Extended Sequence Design (ESD) for tiled, directional pairs. Both derive strand-specific restriction products and block release unless in-silico annealing and ordered ligation reconstruct both intended strands exactly. A staged Design→Hybridization→Restriction-cloning workflow exposes unpaired bases as cohesive ends in compact and nucleotide-level double-strand views, verifies complementarity with a digested vector and requires explicit ligation before product analysis. The same record supports peptide back-translation with explicit initiator handling, editable annotation, PCR/digest/gel simulation and quality-aware Sanger validation. The separated engine, API and interface passed 1,113, 265 and 91 tests, respectively. In an experimentally grounded case study, SSD regenerated insulin glargine A- and B-chain constructs; all four oligonucleotides exactly matched archived synthesis records. NdeI/XhoI cloning into pET-21a(+) preserved the intended open reading frames and generated 5,490-bp and 5,523-bp recombinant plasmids. G-Synth assembled the four author-designated chromatograms into 100%-covered, 100%-identical consensus sequences. Bidirectional overlap was {a_consensus['bidirectional_overlap_percent']:.1f}% for A and {b_consensus['bidirectional_overlap_percent']:.1f}% for B, with 100% F/R agreement; independent Biopython analysis reproduced complete coverage and identity. The exact SSD/ESD release-gate combination was not identified in prior design-automation literature and is, to our knowledge, first automated and disclosed here. Protein expression, folding and bioactivity were not tested."
    )
    add_text(doc, "Keywords: synthetic biology; gene synthesis; restriction cloning; Sanger sequencing; software validation; insulin glargine.")

    add_figure(doc, FIG / "Figure_1_GSynth_workflow_and_architecture.png",
        "Figure 1. G-Synth workflow and software architecture. The design, build, clone and verify stages share a deterministic Python engine. The API stores project provenance and the web workspace renders the engine outputs without reimplementing biological calculations.")

    doc.add_heading("1. Introduction", level=1)
    add_text(doc,
        "Synthetic-biology design–build–test–learn cycles depend on an unbroken relationship between a digital construct and the molecules handled at the bench [13,58–60]. GeneDesign automates synthetic-gene design [61,62], j5 optimizes DNA-assembly plans [58], basicsynbio couples BASIC constructs to transferable build instructions [59], and DNA-BOT links a standardized scheme to laboratory automation [60]. Geneious provides broad commercial cloning and sequence analysis [65], whereas DIVA, CloneCoordinate and InSillyClo extend design–implementation–verification or campaign-level construction planning [69–71]. These systems establish that G-Synth cannot credibly be called the first general end-to-end DNA-construction platform. A distinct gap nevertheless remains for laboratories that order short coding constructs as complementary oligonucleotides: restriction sites may be treated as literal tails instead of strand-specific cut geometry; an enzyme-retained start codon may be duplicated; a reverse molecule may receive the wrong cohesive-end remainder; or an apparently plausible chromatogram may be reported as complete verification despite incomplete evidential support."
    )
    add_text(doc,
        "These problems require a workflow in which each exported sequence is derived from explicit molecular rules and then re-read from the simulated duplex. G-Synth therefore defines SSD for constructs released as one forward/reverse synthesis pair and ESD for longer constructs decomposed into multiple hybridized pairs joined by unique, directional internal overhangs. Their shared end-to-end unit is a traceable molecular record: biological target → synthesis-ready order molecules → reconstructed duplex → recombinant plasmid → annotated construct → sequencing evidence. The engine derives enzyme products from recognition and cut positions, constructs the full coding cassette, creates every strand, re-anneals and re-ligates the strands in silico, checks the outer ends and coding frame, simulates vector digestion and ligation, and maps capillary reads back to the same reference."
    )
    add_text(doc,
        "Insulin glargine A- and B-chain constructs offered a stringent retrospective case study because their long complementary oligonucleotides had been independently synthesized, hybridized, cloned into pET-21a(+), transformed into Escherichia coli DH5α and sequenced. The original draft used Benchling for reverse-complement operations and Geneious Prime for visual cloning and alignment. Here, those inputs and primary trace files were reanalysed with G-Synth as the primary system. The objectives were to (i) reproduce the exact archived orderable molecules from the optimized chain sequences; (ii) validate terminal-end geometry, junctions, open reading frames and plasmid context; (iii) assemble the physical F/R trace evidence and distinguish complete consensus coverage from bidirectional overlap; and (iv) compare G-Synth's focused, auditable workflow with established software categories."
    )

    doc.add_heading("2. Results: software design and implementation", level=1)
    doc.add_heading("2.1 Separation of scientific logic, service layer and interface", level=2)
    add_text(doc,
        "G-Synth comprises three layers (Figure 1). The gsynth_engine Python package contains sequence normalization, codon optimization, nearest-neighbour thermodynamics, enzyme definitions, duplex construction, multi-oligonucleotide assembly, PCR, cloning, ligation, annotation, gel simulation and read verification. It has no web-framework dependency. A Django REST API validates requests, applies access control and serializes immutable calculation outputs together with project provenance. A React/TypeScript workspace renders the outputs, supports editable annotations and guides the user through Design, PCR, Clone, Check, Compare and Learn. This separation makes the engine importable in scripts and permits molecular regression tests to execute without a database or browser."
    )
    add_text(doc,
        "Codon optimization is explicitly conditioned on the selected expression host. Fifteen profiles cover bacterial (Escherichia coli, Bacillus subtilis, Pseudomonas putida, Lactococcus lactis, Corynebacterium glutamicum and Streptomyces coelicolor), yeast (Saccharomyces cerevisiae, Komagataella phaffii, Kluyveromyces lactis and Yarrowia lipolytica), mammalian (human and Chinese-hamster species proxies), insect (Spodoptera frugiperda and Drosophila melanogaster species proxies) and plant (Nicotiana benthamiana) systems. Each profile preserves all 64 raw codon counts, NCBI taxon, source dataset, sample size, GC and snapshot checksum. Counts from the September 2021 FDA HIVE-CUTs/CoCoPUTs release are normalized within each synonymous family [83,84]; RefSeq species aggregates are preferred, with a disclosed GenBank fallback when that snapshot lacks a RefSeq table. Every candidate is translated after optimization, and amino-acid identity is a blocking invariant. A species-wide result is labelled profile-relative CAI rather than strict CAI or predicted yield; a documented highly expressed reference-gene set can override it for strain-, tissue- or cell-line-specific analysis [82,85–87]."
    )
    add_text(doc,
        "The same workflow accepts a peptide and performs host-conditioned reverse translation. Its start logic distinguishes the coding ORF from a mature or post-translationally processed product [88–90]. Automatic mode treats an N-terminal methionine as an existing initiator and otherwise preserves the supplied peptide as a mature product whose initiation is provided upstream by Design. Because sequence alone cannot establish biological maturity, the user may override this inference. Declaring a non-M peptide as a complete expression protein adds exactly one initiator methionine/ATG; declaring it mature never changes its residue sequence. The engine records the input peptide, resolved role and any added initiator, re-translates the generated DNA and transfers the corresponding coding-state flag to Design. It does not predict in-vivo methionine excision."
    )
    doc.add_heading("2.2 Restriction enzymes as cut geometry", level=2)
    add_text(doc,
        "An enzyme is represented by its recognition sequence and top- and bottom-strand cut offsets, not by a prewritten tail. From these offsets the engine derives overhang sequence, length and 5′/3′ polarity for either end of a linear fragment. This matters because the bases that remain on a forward oligonucleotide depend on whether an enzyme is placed at the left or right end and on which strand is ordered. The selectable table contains 109 non-redundant in-site cut geometries representing 289 commercially available enzyme names from REBASE; isoschizomers remain searchable aliases rather than duplicate map sites [75]. Every geometry is exercised in both terminal positions. HindIII, for example, is rendered from A|AGCTT / TTCGA|A cut coordinates and therefore appears with the correct 5′ cohesive end rather than being lost by a recognition-site-only search. Enzymes with ambiguous sites, out-of-site cleavage or two cuts are outside this model and are not claimed as supported. The enzyme-table checksum is stored in provenance records."
    )
    doc.add_heading("2.3 Coding cassette and PCR-free oligonucleotide design", level=2)
    add_text(doc,
        "In Small Sequence Design (SSD), the default expression cassette is ATG (unless the retained enzyme product supplies it), a flexible linker, 6×His tag, a second linker, a thrombin-recognition sequence and the user insert. A stop codon can be retained or added according to the intended downstream fusion. G-Synth checks translation, internal selected-enzyme sites, sequence composition, repeats and secondary-structure warnings. Melting temperatures use nearest-neighbour thermodynamics under recorded ionic and strand-concentration assumptions [80], not a GC-only formula."
    )
    add_text(doc,
        "Within the SSD single-pair length limit, G-Synth generates one forward and one reverse strand with the exact complementary core and asymmetric terminal cohesive ends. Extended Sequence Design (ESD) divides longer constructs into forward/reverse hybridized fragments joined through orthogonal 4–8-nt internal junctions. Every SSD or ESD proposal is simulated end to end: each pair is annealed; ESD fragments are ligated in the declared order; both reconstructed strands must equal the target; every internal overhang must be unique against its reverse complement; and the two outer products must equal the selected restriction-enzyme ends. Export remains blocked if any invariant fails."
    )
    doc.add_heading("2.4 Hybridization as the gate between design and cloning", level=2)
    add_text(doc,
        "Alignment and hybridization answer different questions and are therefore presented as separate modes in one workspace. Alignment may introduce gaps to compare sequence similarity; hybridization instead tests an ungapped antiparallel placement of two molecules entered 5′→3′. G-Synth reverse-complements the second input only for the physical drawing, marks every complementary pair and mismatch, and leaves terminal displacement visible as a strand-specific 5′ or 3′ overhang. The simple view summarizes the duplex core and both ends, whereas the detailed view shows every base, coordinate and polarity. Melting temperature is reported only for a perfectly matched overlap under explicit strand, Na⁺, Mg²⁺ and temperature conditions; it is withheld for mismatched duplexes because the implemented nearest-neighbour model does not include mismatch parameters [80]."
    )
    add_text(doc,
        "A verified SSD or ESD design can be transferred directly to this view with its left and right enzyme identities. Editing either enzyme invalidates the inherited molecules and requires redesign, preventing a sequence from being relabelled with ends it does not possess. Only an exactly complementary duplex can continue to cloning. This implements a visible state transition from an orderable pair, through simulated annealing, to an insert whose exposed ends can be compared with the independently digested vector."
    )
    add_figure(doc,
        INTERFACE / "Hybridization_Workflow.png",
        "Figure 2. G-Synth hybridization view after direct transfer from SSD. Both order molecules remain entered 5′→3′; the lower strand is drawn antiparallel. The 123-bp core is completely paired, while the NdeI-derived 5′-TA and XhoI-derived 5′-TCGA products remain visibly unpaired as cohesive ends.",
    )

    doc.add_heading("2.5 Cloning, annotation and virtual experiments", level=2)
    add_text(doc,
        "The cloning module imports FASTA, GenBank and SnapGene DNA files and includes verified pET-21a(+) and pET-21(+) sequences. It checks vector identity, requires each selected enzyme to cut once and presents the workflow as Vector + duplex Insert → Product. Before ligation, six compatibility checks report strand pairing, unique vector cuts, orientation, site regeneration and reading-frame consequences. Simple and detailed junction views place each exposed insert end beside the complementary vector end. The product map, diagnostic gel and exports remain hidden until the user explicitly selects Ligate compatible ends; this button simulates joining the phosphodiester backbone and is not represented as an experimental ligation. The coordinate-level annotated view handles circular-origin crossings, strand direction and codon-aligned translation. Curated exact matches to common promoters, operators, tags, linkers and cleavage motifs are suggestions, not silently accepted claims; users can create, name, move or delete features and preserve edits in GenBank export."
    )
    add_text(doc,
        "PCR simulation separates a primer's annealing 3′ segment from its intentionally unpaired 5′ extension. The extension becomes part of the amplicon only after polymerase copying, a distinction shown explicitly in the hybridization view. Diagnostic digests produce fragment sizes from the simulated molecule, and virtual gels position those products relative to explicit 100-bp, 1-kb and broad-range marker lists. These gels are labelled as in-silico predictions because migration, topology, partial digestion, staining and band intensity are experimental properties."
    )
    add_figure(doc,
        INTERFACE / "Cloning_PreLigation.jpg",
        "Figure 3A. Pre-ligation compatibility gate. G-Synth shows the NdeI and XhoI vector–insert junctions, their opposite polarities and all six passed checks while withholding recombinant-product analyses.",
    )
    add_figure(doc,
        INTERFACE / "Cloning_Ligated_Product.jpg",
        "Figure 3B. Post-ligation state for the glargine A-chain construct. The 5,490-bp product is created only after explicit in-silico ligation, after which the plasmid map, diagnostic digest, gel and exports become available.",
    )

    doc.add_heading("2.6 Quality-aware capillary-trace verification", level=2)
    add_text(doc,
        "G-Synth detects ABI/ABIF and Standard Chromatogram Format (SCF) content, reads called bases, qualities, peak positions and the four fluorescence channels, evaluates both orientations and aligns reads to a linear or circular reference. Forward and reverse calls are transformed into reference orientation and assembled position by position; equal-support conflicts remain N rather than being resolved in favor of the reference. The interface reports assembled-consensus coverage, consensus identity, the fraction supported by both orientations and F/R agreement within that overlap as distinct quantities. Reference-aligned chromatograms display consensus, strand, base qualities, four-channel peaks and mismatches. Optional user-selected quality trimming is retained for prospective workflows, but quality thresholds must be predefined for the relevant instrument, basecaller and study rather than imposed retrospectively. Phred-style quality values are treated as estimates of base-call error probability [78,79]."
    )
    doc.add_heading("2.7 Originality and application scope", level=2)
    add_text(doc,
        "The originality of G-Synth is the formalization and automation of SSD and ESD, not the invention of restriction cloning, codon optimization, oligonucleotide synthesis or Sanger sequencing individually. GeneDesign, j5, DNA Chisel, basicsynbio, DIVA, CloneCoordinate and InSillyClo provide important but differently scoped automation [58–63,69–71]. In the literature reviewed through 2 September 2026, no system was found that defines SSD and ESD as synthesis-order workflows and requires exact reconstruction of both intended strands from their exported molecules before release, while carrying the same hashed design through restriction cloning and F/R capillary consensus. We therefore describe G-Synth as the first disclosed automation of these specifically defined workflows, to our knowledge, rather than as the first general DNA-design or end-to-end construction platform."
    )
    add_text(doc,
        "This focus supports applications in which small or modular DNA constructs must move rapidly from concept to synthesis: therapeutic-peptide precursors; affinity-tag and protease-cleavage cassettes; vaccine or diagnostic antigens; peptide standards and positive controls; recombinant growth factors, cytokine fragments and antimicrobial peptides; enzyme domains; regulatory elements; reporter constructs; and modular synthetic-biology parts. The glargine A/B case study is especially relevant because therapeutic peptides impose simultaneous constraints on exact amino-acid identity, purification architecture, cleavage logic, reading frame and post-cloning verification. G-Synth addresses the DNA design-to-sequence portion end to end. Expression yield, processing, oxidative folding, higher-order structure, potency, safety and regulatory comparability remain separate experimental stages."
    )

    doc.add_heading("3. Validation results", level=1)
    add_text(doc,
        "Validation combined unit, property, API and interface tests with a retrospective experimental case study. The engine suite contains golden examples, error-path tests, host-profile invariants, peptide-start decisions, enzyme-wide properties and antiparallel-hybridization cases covering 5′/3′ overhangs, blunt ends, ambiguity, mismatches and thermodynamic reporting boundaries. The defining assembly property is equality between the intended construct and both strands reconstructed from the exported oligonucleotides. API tests cover request validation, host selection, peptide-to-DNA serialization, authentication, projects, sequence import, hybridization, PCR and security. Interface tests cover the design and validation views, host-profile and peptide-role selection, automatic Design-to-Hybridization-to-Restriction-cloning handoff, simultaneous compact and nucleotide-level duplex evidence, editable annotations, primer-tail hybridization, gels, Learn content and responsive behaviour. The complete suites were executed from the manuscript working tree on 2 September 2026: 1,113/1,113 engine, 265/265 API and 91/91 interface tests passed."
    )
    add_table(doc,
        ["Validation layer", "Primary question", "Result"],
        [
            ["Molecular engine", "Do algorithms preserve molecular invariants and reject invalid inputs?", "1,113 passed"],
            ["HTTP/API", "Are calculations exposed consistently with accounts, projects and security controls?", "265 passed"],
            ["Web interface", "Are outputs rendered and interactions retained without reimplementing biology?", "91 passed"],
            ["Insulin glargine design", "Do archived A/B inputs regenerate the four molecules synthesized at the bench?", "4/4 exact"],
            ["Cloning simulation", "Are NdeI/XhoI ends, junctions and ORFs coherent in pET-21a(+)?", "A and B passed"],
            ["Sequencing evidence", "Do oriented F/R reads assemble across each insert and agree in their overlap?", "100% consensus coverage and identity; 100% overlap agreement"],
        ], [1.25, 3.85, 1.5])

    doc.add_heading("4. Insulin glargine A/B-chain case study", level=1)
    doc.add_heading("4.1 Input sequences and cassette logic", level=2)
    add_text(doc,
        "The optimized A-chain coding insert was 66 bp and encoded GIVEQCCTSICSLYQLENYCG; the B-chain insert was 99 bp and encoded FVNQHLCGSHLVEALYLVCGERGFFYTPKTRR. The second sequence includes the two C-terminal arginine residues characteristic of insulin glargine's B-chain precursor context. For each construct, G-Synth placed the insert after MGSSHHHHHHSSGLVPRGS, retained TAA termination and selected NdeI at the 5′ end and XhoI at the 3′ end. NdeI's retained bases supply the initiating ATG, so adding a separate ATG would duplicate the start. The resulting products were MGSSHHHHHHSSGLVPRGSGIVEQCCTSICSLYQLENYCG and MGSSHHHHHHSSGLVPRGSFVNQHLCGSHLVEALYLVCGERGFFYTPKTRR."
    )
    add_figure(doc, INTERFACE / "Design_Release_Gate.png",
        "Figure 4. G-Synth design view for the glargine A-chain input. The release gate reports exact two-strand reconstruction, terminal-end compatibility and the NdeI start-codon note before displaying the orderable construct.")

    a = design["chains"]["A"]
    b = design["chains"]["B"]
    add_table(doc,
        ["Metric", "A-chain construct", "B-chain construct"],
        [
            ["Optimized chain insert", f"{a['input']['length_bp']} bp", f"{b['input']['length_bp']} bp"],
            ["Forward / reverse oligo", f"{a['design']['forward_length_nt']} / {a['design']['reverse_length_nt']} nt", f"{b['design']['forward_length_nt']} / {b['design']['reverse_length_nt']} nt"],
            ["Forward / reverse GC", f"{a['design']['forward_gc_percent']:.1f} / {a['design']['reverse_gc_percent']:.1f}%", f"{b['design']['forward_gc_percent']:.1f} / {b['design']['reverse_gc_percent']:.1f}%"],
            ["Forward / reverse Tm", f"{a['design']['forward_tm_c']:.1f} / {a['design']['reverse_tm_c']:.1f} °C", f"{b['design']['forward_tm_c']:.1f} / {b['design']['reverse_tm_c']:.1f} °C"],
            ["Observed terminal products", "5′ TA / 5′ TCGA", "5′ TA / 5′ TCGA"],
            ["Exact match to archived synthesis record", "Forward and reverse", "Forward and reverse"],
            ["Recombinant pET-21a(+) length", f"{a['cloning']['recombinant_length_bp']:,} bp", f"{b['cloning']['recombinant_length_bp']:,} bp"],
            ["Junction restriction sites", "NdeI and XhoI regenerated", "NdeI and XhoI regenerated"],
        ], [2.25, 2.18, 2.17])
    add_text(doc,
        "The reverse oligonucleotide is two nucleotides longer than the forward molecule in each pair (127 versus 125 nt for A; 160 versus 158 nt for B). This is an expected consequence of the strand-specific cohesive-end remainders, not a discrepancy. All four generated sequences were identical to the corresponding archived order records, and the complementary cores reconstructed without mismatch."
    )

    doc.add_heading("4.2 pET-21a(+) cloning and feature-level validation", level=2)
    add_text(doc,
        "G-Synth digested the bundled pET-21a(+) sequence with NdeI and XhoI, checked the insert's physical ends and simulated directional ligation. Both junctions were compatible and regenerated their recognition sites. The A and B recombinant molecules were 5,490 bp and 5,523 bp. Translation from the NdeI-supplied ATG reached the intended TAA without frameshift. Because the inserts terminate before the vector's downstream coding region, the pET-21a(+) C-terminal His tag is not translated; the designed N-terminal 6×His cassette remains present. Similarly, cloning at NdeI replaces the vector's optional N-terminal T7 tag with the designed cassette. These are sequence consequences, not predictions of expression, solubility, cleavage efficiency or biological activity."
    )
    add_figure(doc, INTERFACE / "Annotated_Glargine_A.jpg",
        "Figure 5A. G-Synth coordinate-level view of the insulin glargine A-chain recombinant pET-21a(+) construct. Editable feature tracks, codon-aligned translation and the regenerated NdeI and XhoI junctions are shown across the circular origin.", width=6.35)
    add_figure(doc, INTERFACE / "Annotated_Glargine_B.jpg",
        "Figure 5B. G-Synth coordinate-level view of the insulin glargine B-chain recombinant pET-21a(+) construct, using the same editable annotation and strand-aware coordinate system.", width=6.35)

    doc.add_heading("4.3 Experimental synthesis and cloning evidence", level=2)
    add_text(doc,
        "The retrospective data originated from an independent wet-laboratory project at the Genomics Technology Platform of the Higher School of Biological Sciences of Oran. Forward and reverse strands were synthesized on a MerMade 4 instrument at 1 µmol scale, purified, and annealed in 2× SSC with 0.1% SDS after 5 min at 95 °C followed by gradual cooling. The source draft reports single-stranded concentrations of 1,688.51–1,709.40 ng/µL. After hybridization, ScanDrop 260/280 ratios were 1.92 for A and 1.88 for B; Bioanalyzer sizes were 129 bp and 173 bp. These apparent sizes exceed the 123-bp and 156-bp duplex references and should be interpreted within the sizing limitations of the assay rather than as base-resolved sequence evidence."
    )
    add_text(doc,
        "pET-21a(+) was double-digested with NdeI and XhoI, gel-purified and ligated to each hybridized insert at a nominal 3:1 insert:vector molar ratio. Ligation products were transformed into E. coli DH5α and selected on ampicillin. Colony PCR and capillary sequencing were performed as described in the source experimental record. The present software paper reuses these physical observations as a case study but does not claim expression, purification, oxidative refolding, receptor activity or therapeutic equivalence."
    )

    doc.add_heading("4.4 Post-sequencing validation by G-Synth", level=2)
    add_text(doc,
        "Validation was restricted to the four author-designated correct files: A Forward Seq.ab1, A Reverse Seq.ab1, B Forward Seq.ab1 and B Reverse Seq.ab1. No other chromatogram was admitted. The reads were analysed against the 123-bp A and 156-bp B cassette references. Although the filenames used the .ab1 extension, binary magic-number inspection identified SCF content; G-Synth therefore parsed the files by content rather than extension. For each chain, the forward and reverse read were independently placed, reverse-complemented when required, clipped to the reference and merged position by position into a reference-guided consensus. Consensus coverage denotes positions with at least one oriented call; bidirectional overlap denotes positions called from both orientations; overlap agreement denotes the fraction of those positions with identical F/R calls. Because the source experiment did not prespecify a trace-quality threshold, consensus coverage, identity and overlap agreement were the primary retrospective endpoints; a fixed-threshold quality analysis was retained only as supplementary sensitivity information."
    )
    add_table(doc,
        ["Chain", "Reference", "Consensus coverage", "Consensus identity", "Consensus differences", "F/R overlap", "F/R agreement", "Independent coverage"],
        [
            ["A", "123 bp", f"{a_consensus['coverage_percent']:.1f}%", f"{a_consensus['identity_percent']:.1f}%", "0", f"{a_consensus['bidirectional_overlap_percent']:.1f}%", f"{a_consensus['bidirectional_overlap_agreement_percent']:.1f}%", "100.0%"],
            ["B", "156 bp", f"{b_consensus['coverage_percent']:.1f}%", f"{b_consensus['identity_percent']:.1f}%", "0", f"{b_consensus['bidirectional_overlap_percent']:.1f}%", f"{b_consensus['bidirectional_overlap_agreement_percent']:.1f}%", "100.0%"],
        ], [0.38, 0.58, 0.82, 0.78, 0.82, 0.67, 0.78, 0.82], font_size=7.3)
    add_text(doc,
        f"The assembled consensus covered 100% of both references with 100% identity and no consensus differences. A contained {a_consensus['bidirectional_overlap_percent']:.1f}% bidirectional overlap and B contained {b_consensus['bidirectional_overlap_percent']:.1f}%; all overlapping F/R calls agreed (100%). Thus the paired reads jointly span each complete insert reference, while the overlap is independently observed on both orientations. Independent local alignments in Biopython {biopython_version} reproduced 100% combined coverage and 100% identity across aligned bases. The archived Geneious Prime figures reported visually concordant alignments but were not re-executed in the scripted environment."
    )
    add_callout(doc,
        "Scientific interpretation",
        "The oriented F/R assembly establishes complete consensus coverage, 100% identity and perfect agreement wherever the two reads overlap. This validates the intended A- and B-chain insert sequences at consensus-call level. It does not imply that every position was independently observed from both orientations, and the insert-only references do not test both vector–insert junctions; a prospective validation study should address those additional endpoints with predefined sequencing criteria.",
        fill=CAUTION,
    )
    add_figure(doc, FIG / "Figure_4_GSynth_reference_aligned_viewer.jpg",
        "Figure 6. G-Synth's interactive reference-aligned chromatogram view for representative approved A- and B-chain reads. Base calls, consensus, quality context and four-channel peaks are displayed in reference orientation.")
    add_figure(doc, EVIDENCE / "sequencing_validation" / "Figure_GSynth_Sanger_Approved_Traces.png",
        "Figure 7. G-Synth F/R consensus validation from the four author-designated chromatograms. The oriented reads jointly cover 100% of both references with 100% consensus identity. Green shading marks bidirectional overlap (56.1% for A and 71.2% for B); every overlapping F/R call agreed.")

    doc.add_heading("5. Relationship to existing software", level=1)
    add_text(doc,
        "The comparison below separates G-Synth's precise contribution from adjacent functionality. GeneDesign and DNA Chisel focus on synthetic-sequence design and optimization [61–63]; j5, basicsynbio and DNA-BOT address assembly planning or build automation [58–60]; Geneious is a broad commercial environment [65]; DIVA, CloneCoordinate and InSillyClo cover increasingly integrated construction workflows [69–71]; and Biopython supplies reusable analysis primitives [66]. G-Synth is narrower than several of these systems. Its contribution is the coupled SSD/ESD synthesis-order model, exact two-strand release gate and continuity into post-Sanger consensus."
    )
    add_table(doc,
        ["Tool/category", "Primary scope", "Relationship to SSD/ESD claim"],
        [
            ["G-Synth", "SSD/ESD order molecules, restriction geometry, exact reconstruction, cloning and F/R consensus", "Implements the claimed workflows and blocks export when molecular invariants fail"],
            ["GeneDesign / DNA Chisel", "Synthetic-gene design and sequence optimization", "Foundational prior art for design automation; no identical SSD/ESD release gate identified"],
            ["j5 / basicsynbio / DNA-BOT", "Assembly design, standardized parts and laboratory build instructions", "Broader or method-specific assembly automation; not the same synthesis-order abstraction"],
            ["Geneious Prime", "Broad commercial cloning, annotation and sequence analysis", "Secondary qualitative comparator; no automated head-to-head benchmark was executed"],
            ["DIVA / CloneCoordinate / InSillyClo", "Integrated construction planning, coordination or Golden Gate/MoClo campaigns", "Directly limits any general end-to-end novelty claim; distinct from the exact SSD/ESD gate"],
            ["Biopython", "Programmatic parsing, alignment and molecular-biology primitives", "Independent executable numerical comparator; equivalent logic can be scripted"],
        ], [1.45, 2.55, 2.6], font_size=7.7)
    add_text(doc,
        "This table is a workflow-positioning comparison, not a speed or accuracy benchmark. No proprietary Geneious algorithm was reverse-engineered, and no claim is made that G-Synth has broader functionality. Instead, the archived Geneious views provide an independent qualitative check, while the numerical cross-check uses openly scriptable Biopython alignment. A future prospective benchmark should use blinded constructs, identical trimming thresholds and predefined acceptance criteria across platforms."
    )

    doc.add_heading("6. Discussion", level=1)
    add_text(doc,
        "The insulin case study demonstrates a useful distinction between design correctness, complete assembled coverage and bidirectional support. The exact match between four G-Synth-generated oligonucleotides and the archived synthesis records, together with correct duplex reconstruction and pET-21a(+) junctions, establishes that the application's design logic reproduces the intended molecular plan. G-Synth's F/R assembly yields complete consensus coverage, 100% identity and 100% agreement across the bidirectional overlaps. Positions outside those overlaps are supported by one oriented read rather than two, a distinction the software reports explicitly instead of obscuring within a single coverage value."
    )
    add_text(doc,
        "Restriction-aware design is particularly vulnerable to plausible-looking errors because recognition sites are visually familiar while cleavage products are strand-specific. Storing cut positions and calculating remainders avoids pair-specific tail constants. Reading the ends back from the reconstructed molecule adds an independent internal check. The same philosophy applies to primer extensions and virtual gels: G-Synth shows what the model establishes and labels what remains experimental."
    )
    add_text(doc,
        "Several limitations define the scope. First, long-oligonucleotide chemical synthesis can produce truncations and heterogeneous products; computational exactness does not predict synthesis yield or purification success. Second, orthogonal internal overhangs reduce one class of misassembly but do not model ligation kinetics. Third, restriction-cloning simulation does not predict methylation sensitivity, star activity, buffer compatibility, partial digestion or topology-dependent migration. Fourth, codon optimization and ORF integrity do not establish expression, folding, disulfide formation or insulin receptor activity. Separate production of A and B chains requires validated purification, cleavage, oxidative assembly, structural characterization and bioactivity studies [6–9,35,36,45–47,50,54–57]. Finally, the retrospective traces require prospective confirmation with predefined quality criteria and outside-vector primers spanning both junctions."
    )
    add_text(doc,
        "The appropriate next validation step is a preregistered, prospective design-to-sequence study. At least three independent constructs should include the NdeI/XhoI glargine path, a mixed 5′/3′ cohesive-end path and an internal-site negative control. Required records should include exported oligos, enzyme and vector lots, uncut/single-cut/double-cut controls, colony screening, raw capillary files, full junction and insert coverage, and independent review. User-centred validation with bench scientists should be reported separately from algorithmic correctness."
    )

    doc.add_heading("7. Conclusions", level=1)
    add_text(doc,
        "G-Synth formalizes Small Sequence Design and Extended Sequence Design as auditable synthesis-order workflows and carries the same digital molecule through peptide back-translation or host-selected codon optimization, oligonucleotide ordering, exact reconstruction, explicit hybridization, staged ligation, editable annotation, virtual PCR/digest/gel experiments and quality-aware post-sequencing review. The software passed 1,469 automated tests across its engine, API and interface. In the insulin glargine case study, SSD exactly regenerated all archived A/B oligonucleotides, produced coherent NdeI/XhoI recombinant pET-21a(+) designs and assembled the approved F/R chromatograms into 100%-covered, 100%-identical consensus sequences with perfect agreement in their bidirectional overlaps. The central contribution is not generic end-to-end software, but the first disclosed automation, to our knowledge, of the specifically defined SSD/ESD two-strand release gate and its continuity into post-sequencing evidence."
    )

    doc.add_heading("8. Materials and Methods", level=1)
    doc.add_heading("8.1 Reproducible software analyses", level=2)
    add_text(doc,
        f"The design and cloning case study was executed with tools/publication/design_glargine_case_study.py. The script records the source-manuscript SHA-256 hash, regenerates the oligonucleotides, asserts exact equality with the archived Table 3 molecules, translates the cassette, simulates pET-21a(+) cloning and writes publication_evidence/glargine_ab_design_and_cloning.json. Trace reanalysis was executed with tools/publication/validate_insulin_correct_traces.py, which admits only the four author-designated files, hashes every reference and trace, constructs oriented G-Synth F/R consensus sequences, quantifies bidirectional overlap and agreement, performs an independent Biopython {biopython_version} local alignment and writes a JSON evidence record plus figure. Full commands are supplied in the Supporting Information."
    )
    add_text(doc,
        "The codon profiles were reconstructed with tools/update_hive_codon_tables.py from FDA HIVE service object 537 and committed as gsynth_engine/data/codon_usage_hive_2021.json. The updater requests genomic species data with descendant taxa, prefers RefSeq and falls back to GenBank only for K. phaffii and N. benthamiana, validates the returned taxon, requires exactly 64 non-negative codon counts and checks that their sum equals the reported total. The application divides each raw count by the maximum count among synonymous codons for that amino acid. The shipped September 2021 snapshot was retrieved on 2 September 2026 and is identified by SHA-256 in the API. Because these are species-wide genomic profiles rather than prespecified highly expressed genes, the geometric-mean score is reported as profile-relative CAI. No expression-yield inference is made. The machine-readable validation record publication_evidence/codon_host_profile_validation.json confirms 15 complete and numerically distinct 64-codon weight vectors and preservation of a 100-residue benchmark protein under every profile."
    )
    add_text(doc,
        "Peptide-start and restriction-catalogue checks were regenerated with tools/publication/validate_peptide_and_enzyme_logic.py. The resulting publication_evidence/peptide_and_enzyme_validation.json records four automatic and overridden peptide-role cases, host-wise reverse translation of the mature glargine A chain, translation hashes, the 109/289 geometry-to-name inventory and the independently derived HindIII 5′-AGCT end. The mature chain remained unchanged and did not receive an artificial ATG under all 15 host profiles; the complete-ORF override added exactly one initiator when required."
    )
    doc.add_heading("8.2 Computational hybridization and staged ligation", level=2)
    add_text(doc,
        "Both hybridization inputs were normalized as unambiguous DNA written 5′→3′. The second was reverse-complemented for offset evaluation, and every ungapped antiparallel placement was scored as twice the number of complementary columns minus three times the number of mismatches. Deterministic ties favoured more pairs, fewer mismatches, longer overlap, fewer terminally unpaired bases and the placement closest to flush. Internal gaps were not introduced because they represent bulges rather than terminal cohesive ends. Unpaired terminal intervals were classified by strand, side, polarity and sequence; overlaps shorter than four complementary bases were reported as insufficient. Perfect overlaps were evaluated with the implemented SantaLucia nearest-neighbour and Owczarzy salt-correction model under recorded buffer conditions; mismatched overlaps were not assigned a melting temperature."
    )
    add_text(doc,
        "For cloning handoff, the inherited duplex, enzyme pair and provenance were transferred without rewriting the strands. The vector was independently cut at the two selected unique sites. G-Synth compared the polarity and reverse-complement sequence of both insert–vector end pairs, then evaluated duplex pairing, cut uniqueness, orientation, regenerated sites and coding frame. Recombinant maps and downstream virtual analyses were enabled only after an explicit in-silico ligation state transition."
    )
    doc.add_heading("8.3 Biological inputs, cassette design and G-Synth analysis", level=2)
    add_text(doc,
        "Insulin glargine A- and B-chain amino-acid sequences were retrieved from DrugBank record DB00047 (database version 5.1.10). The primary experimental study used GenScript reverse translation and the VectorBuilder codon-optimization tool for Escherichia coli; the exact optimized nucleotide inputs were preserved unchanged. Each cassette comprised an initiating ATG, the N-terminal peptide MGSSHHHHHHSSGLVPRGS (6×His tag, linkers and thrombin-recognition segment), the optimized glargine chain and TAA. G-Synth SSD took the preserved optimized sequence as input, derived NdeI- and XhoI-compatible strand products from cut geometry, generated the order molecules, re-annealed them in silico, verified exact reconstruction on both strands, translated the reconstructed coding sequence, and simulated directional insertion into pET-21a(+). Benchling and ExPASy operations reported in the primary draft are retained as historical design provenance; the present validation was performed first in G-Synth."
    )
    doc.add_heading("8.4 Oligonucleotide synthesis, deprotection and purification", level=2)
    add_text(doc,
        "Forward and reverse A- and B-chain oligonucleotides were synthesized at the ESSBO Genomics Technology Platform by phosphoramidite chemistry on a MerMade 4 synthesizer (LGC Biosearch Technologies, Hoddesdon, UK) using Universal Support Columns, DMT-off (N-iPr), 1000 Å and 1 µmol scale. The source program comprised acetonitrile initialization washes; two deblocking deliveries of 3% trichloroacetic acid in dichloromethane; activator/amidite coupling in two subinjection stages; acetonitrile washing; Cap A and Cap B deliveries; iodine oxidation; final deblocking and washing; and cleavage/deprotection in 1 mL concentrated ammonia for 18 h at 70 °C. All volumes, durations and equalization times are reproduced without omission in Table S5."
    )
    add_text(doc,
        "Purification used an equal volume of butanol prechilled to −20 °C, 30 s vortexing, 30 min at −20 °C and centrifugation at 14,000 rpm for 15 min at 4 °C. Pellets were washed twice with 1 mL of 70% ethanol at −20 °C; each wash was followed by 10 min at 14,000 rpm and 4 °C. Pellets were air-dried for 10–15 min and resuspended in ultrapure sterile water. Purity and concentration were assessed by ScanDrop spectrophotometry (Analytik Jena) and 2% agarose gel electrophoresis at 100 V for 45 min. Rotor radius was not recorded; consequently rpm cannot be converted reproducibly to relative centrifugal force."
    )
    doc.add_heading("8.5 Hybridization and analytical quality control", level=2)
    add_text(doc,
        "The 20× SSC stock contained 3 M NaCl and 0.3 M sodium citrate and was adjusted to pH 7.0. For each duplex, 10 µL of each purified complementary strand was combined with 20 µL of 2× SSC containing 0.1% SDS, heated at 95 °C for 5 min and cooled gradually to room temperature over 30 min. Samples were stored at 4 °C for short-term use or −80 °C for long-term storage. Duplex concentration, purity and apparent size were assessed by ScanDrop, 2% agarose electrophoresis (100 V, 45 min) and an Agilent 2100 Bioanalyzer with DNA 1000 reagents (25–1000-bp stated range; ±5% sizing resolution; 15% quantitation CV for 100–500 bp)."
    )
    doc.add_heading("8.6 Vector preparation, ligation, transformation and clone screening", level=2)
    add_text(doc,
        "pET-21a(+) was digested in 50 µL containing 1 µg vector, 1 µL NdeI (10 U/µL), 1 µL XhoI (10 U/µL), 5 µL 10× CutSmart buffer and nuclease-free water for 2 h at 37 °C. The digest was resolved on 1% agarose and gel-purified with an Invitrogen kit. Ligation used a 3:1 insert:vector molar ratio calculated with NEBioCalculator, 50 ng digested vector, the calculated insert quantity, 1 µL T4 DNA ligase (400 U/µL), 2 µL 10× ligase buffer and water to 20 µL; reactions were held at 16 °C for 2 h and 4 °C overnight. Ligation products were assessed by ScanDrop, Bioanalyzer and 2% agarose electrophoresis (100 V, 45 min)."
    )
    add_text(doc,
        "Competent E. coli DH5α cells were heat-shocked for 50 s at 42 °C, recovered in SOC medium and plated on LB agar containing 100 µg/mL ampicillin. Colonies were screened by PCR using both T7 and insert-specific primer pairs. Positive clones were subjected to plasmid miniprep (Thermo Fisher Scientific) and capillary sequencing. LB medium and ampicillin were from Sigma-Aldrich; SOC medium and DH5α cells were from Thermo Fisher Scientific."
    )
    doc.add_heading("8.7 Capillary sequencing, G-Synth consensus and secondary confirmation", level=2)
    add_text(doc,
        "Each 20-µL BigDye Terminator v3.1 reaction contained 20 ng purified DNA, 2 µL BigDye Terminator v3.1, 3 µL 5× sequencing buffer and 1 µL forward or reverse primer at 3 µM. Cycling comprised 96 °C for 60 s followed by 35 cycles of 96 °C for 10 s, 50 °C for 5 s and 60 °C for 4 min. Products were ethanol-precipitated and analysed on an Applied Biosystems 3500 Genetic Analyzer with POP-7 polymer at the ESSBO Genomics Technology Platform."
    )
    add_text(doc,
        "Only A Forward Seq.ab1, A Reverse Seq.ab1, B Forward Seq.ab1 and B Reverse Seq.ab1 were admitted to the present analysis. Content signatures identified the files as SCF despite their extensions. G-Synth parsed calls, qualities, peak locations and four-channel traces; evaluated both orientations; transformed placed reads into reference orientation; and merged them position by position. Consensus calls were selected by supporting-read count and then cumulative quality; exact ties were emitted as N. Consensus coverage, identity, bidirectional overlap and overlap agreement were computed separately. Because the source experiment had no prespecified trace-quality acceptance threshold, no post hoc cutoff was used for the primary retrospective endpoint. Biopython PairwiseAligner supplied an independent executable numerical comparison, and the archived Geneious Prime 2024.0.7 figure using these same four approved files supplied a secondary qualitative confirmation."
    )

    doc.add_heading("Data and Software Availability", level=1)
    add_text(doc,
        "The G-Synth source code, tests and documentation are available under the MIT License at https://github.com/Midotech31/g-synth-app. A versioned archival DOI will be inserted after repository release and Zenodo deposition. Reproducibility scripts, the immutable primary-draft record and JSON evidence are under tools/publication/ and publication_evidence/. Glargine records PQ362993 and PQ362994 are in GenBank. Raw capillary files contain no human data; deposition details will be finalized before submission or supplied confidentially to reviewers when permitted."
    )
    doc.add_heading("Supporting Information", level=1)
    add_text(doc,
        "Complete designs and oligonucleotides; algorithmic invariants and executable commands; experimental and QC records; trace manifest; interface evidence; comparison notes; correction ledger; and release checklist (DOCX and PDF)."
    )
    doc.add_heading("Author Contributions", level=1)
    add_text(doc,
        "M.M. conceived G-Synth and the PCR-free construct strategy, defined the scientific requirements, supervised software development and validation, supplied and adjudicated the experimental evidence, and drafted the manuscript. Z.Y.Z. contributed experimental investigation and manuscript review. Y.S., A.I.H., M.A., K.B., H.S.B. and S.B. contributed experimental investigation. D.S. contributed supervision, experimental review and manuscript review."
    )
    doc.add_heading("Funding and Acknowledgments", level=1)
    add_text(doc,
        "This work was conducted at the Genomics Technology Platform of the Higher School of Biological Sciences of Oran and supported by the General Directorate for Scientific Research and Technological Development (DGRSDT), Ministry of Higher Education and Scientific Research of Algeria."
    )
    doc.add_heading("Generative-AI use", level=2)
    add_text(doc,
        "OpenAI Codex (GPT-5 family, Codex desktop application; accessed June–September 2026) assisted under Mohamed Merzoug's direct supervision with software review and refactoring, test scaffolding and execution, documentation and manuscript drafting, literature-discovery assistance, copy-editing, release preparation and document-layout quality control. The authors supplied and interpreted the experimental evidence, adjudicated the biological assumptions, reviewed and edited the AI-assisted output and retain responsibility for the accuracy, originality, licensing and scientific claims. Generative AI was not treated as an author, did not generate experimental observations and did not independently determine novelty or scientific validity."
    )
    doc.add_heading("Conflict of Interest", level=1)
    add_text(doc, "The authors declare no competing financial interest.")

    doc.add_heading("References", level=1)
    add_reference_list(doc, main_references(primary))
    PUBLICATION.mkdir(parents=True, exist_ok=True)
    doc.save(MAIN_OUT)


def build_si(design: dict, seq: dict, primary: dict) -> None:
    doc = Document()
    biopython_version = seq.get("analysis", {}).get("biopython_version", "version recorded in the evidence JSON")
    configure_document(doc, "G-Synth — Supporting Information")
    add_title_block(doc, "Supporting Information",
        "G-Synth: Small Sequence Design and Extended Sequence Design for auditable synthesis-ready DNA and post-sequencing validation")
    add_authors(doc)
    add_text(doc,
        "Contents: complete glargine sequences; SSD/ESD algorithmic invariants; reproducibility commands; unabridged experimental methods, synthesis script and QC record from the primary ACS draft; trace-file manifest; G-Synth evidence; secondary Geneious comparisons; and a scientific correction ledger."
    )

    doc.add_heading("S1. Exact glargine design inputs and outputs", level=1)
    add_table(doc,
        ["Sequence type", "Description", "Sequence (preserved from primary draft)"],
        primary["tables"][2][1:], [1.15, 1.55, 3.9], font_size=7.7)
    for chain in ("A", "B"):
        d = design["chains"][chain]
        doc.add_heading(f"S1.{1 if chain == 'A' else 2} Insulin glargine {chain}-chain", level=2)
        add_text(doc, f"Optimized coding input ({d['input']['length_bp']} bp):")
        add_callout(doc, "5′→3′ input", d["input"]["optimised_coding_sequence"], fill=LIGHT)
        add_text(doc, f"Expected chain peptide: {d['input']['expected_chain_peptide']}")
        add_text(doc, f"Forward order molecule ({d['design']['forward_length_nt']} nt):")
        add_callout(doc, "5′→3′ forward", d["design"]["forward_oligo_5_to_3"], fill=PALE_TEAL)
        add_text(doc, f"Reverse order molecule ({d['design']['reverse_length_nt']} nt):")
        add_callout(doc, "5′→3′ reverse", d["design"]["reverse_oligo_5_to_3"], fill=PALE_GREEN)
        add_text(doc, f"Translated cassette: {d['design']['translated_product']}")
        add_text(doc, "Validation flags: exact archived forward match; exact archived reverse match; exact duplex core; expected product translation.")

    add_table(doc,
        ["Order molecule", "Nucleotide sequence (5′→3′)", "Accession / translated product"],
        primary["tables"][3][1:], [1.35, 3.75, 1.5], font_size=6.8)

    doc.add_heading("S2. Molecular invariants enforced by SSD and ESD", level=1)
    add_bullets(doc, [
        "Input sequences are normalized to unambiguous DNA and coding designs are translated before release.",
        "SSD releases one complementary forward/reverse synthesis pair; ESD tiles longer targets into multiple ordered complementary pairs with directional internal junctions.",
        "Recognition sequences and strand cut offsets determine terminal products; overhangs are not stored as pair-specific oligo literals.",
        "An enzyme supplies an initiating ATG only when the bases retained after cleavage place ATG at the coding boundary.",
        "The forward and reverse oligonucleotide cores must be exact reverse complements.",
        "Annealing and ordered ligation must reconstruct the designed molecule on both strands, base for base.",
        "Internal assembly overhangs must be unique against all used overhangs and their reverse complements; designs widen beyond four bases when required.",
        "Outer reconstructed ends must equal the selected enzyme products with correct sequence and polarity.",
        "A cloning vector must contain exactly one site for each selected enzyme; physical vector and insert ends must anneal.",
        "Coding-junction translation, stop-codon consequences, vector-tag outcomes and site regeneration are reported explicitly.",
        "Full sequence verification requires complete requested-region coverage and zero confident unexplained differences.",
    ])

    doc.add_heading("S3. Reproducibility", level=1)
    add_text(doc, "Execute from the repository root with Python 3.12.13 and the locked project dependencies:")
    commands = (
        "python tools/publication/design_glargine_case_study.py --manuscript <source-draft.docx> --output publication_evidence/glargine_ab_design_and_cloning.json\n"
        "python tools/publication/validate_insulin_correct_traces.py --root <reference-data-root> --trace-dir publication_evidence/validated_traces --output-dir publication_evidence/sequencing_validation\n"
        "python -m pytest gsynth_engine/tests -q\n"
        "cd django_app && python -m pytest -q\n"
        "cd ../frontend && npm test -- --run"
    )
    add_callout(doc, "Commands", commands, fill=LIGHT)
    add_table(doc,
        ["Evidence object", "Location", "Purpose"],
        [
            ["Design/cloning JSON", "publication_evidence/glargine_ab_design_and_cloning.json", "Exact molecules, segments, translations, junctions, warnings and source hash"],
            ["Sequencing JSON", "publication_evidence/sequencing_validation/glargine_approved_trace_validation.json", "Approved-file policy, hashes, primary F/R consensus metrics, post hoc quality sensitivity and Biopython check"],
            ["Design script", "tools/publication/design_glargine_case_study.py", "Executable glargine case study with equality assertions"],
            ["Trace script", "tools/publication/validate_insulin_correct_traces.py", "Executable four-file sequencing analysis and figure generation"],
        ], [1.25, 2.8, 2.55], font_size=8.2)

    doc.add_heading("S4. Unabridged experimental oligonucleotide synthesis and cloning record", level=1)
    add_text(doc,
        f"The record below was extracted from the primary ACS draft {primary['source_filename']} (SHA-256 {primary['source_sha256']}) without removing technical steps. Laboratory notebooks and instrument exports remain the authoritative operational record. Oligonucleotides were synthesized on a MerMade 4 instrument using 1-µmol, 1000-Å universal supports and DMT-off processing."
    )
    add_text(doc,
        "Purification used an equal volume of cold butanol, 30 s vortexing, 30 min at −20 °C and 15 min centrifugation at 14,000 rpm and 4 °C. Pellets were washed twice with 1 mL 70% ethanol at −20 °C, centrifuged for 10 min, air-dried and resuspended. Exact rotor radius is not recorded, so relative centrifugal force cannot be reconstructed from rpm; future protocols should report ×g."
    )
    add_text(doc,
        "Hybridization combined 10 µL of each complementary strand with 20 µL 2× SSC containing 0.1% SDS. Samples were heated at 95 °C for 5 min and cooled to room temperature over 30 min. pET-21a(+) digestion used 1 µg vector, 10 U NdeI, 10 U XhoI and 1× CutSmart buffer in 50 µL for 2 h at 37 °C. Gel-purified vector was ligated to duplex insert at a nominal 3:1 molar ratio with 50 ng vector and 400 U T4 DNA ligase in 20 µL at 16 °C for 2 h followed by 4 °C overnight. E. coli DH5α cells were heat-shock transformed for 50 s at 42 °C, recovered in SOC and selected on 100 µg/mL ampicillin."
    )
    add_table(doc,
        ["Measurement", "A forward", "A reverse", "B forward", "B reverse"],
        [
            ["ssDNA concentration (ng/µL)", "1702.52", "1709.40", "1688.51", "1703.60"],
            ["dsDNA concentration by ScanDrop (ng/µL)", "830", "830", "1936", "1936"],
            ["dsDNA A260/A280", "1.92", "1.92", "1.88", "1.88"],
            ["Bioanalyzer apparent size (bp)", "129", "129", "173", "173"],
            ["Bioanalyzer concentration (ng/µL)", "1.15", "1.15", "3.39", "3.39"],
        ], [1.85, 1.18, 1.18, 1.18, 1.18], font_size=8.2)
    add_text(doc,
        "Caution: UV absorbance, gel mobility and Bioanalyzer sizing assess concentration, purity and approximate fragment size. They do not establish base identity. The different concentration values from ScanDrop and Bioanalyzer reflect method and dilution context and should not be combined without the original dilution records."
    )
    add_figure(doc, SOURCE_FIG / "Figure_S2_WetLab_Oligonucleotide_QC.png",
        "Figure S2. Primary wet-laboratory quality-control evidence. Agarose electrophoresis of the synthesized strands and hybridized duplexes, and Bioanalyzer electropherograms, are reproduced from the source experimental record. These assays support presence and approximate size, not nucleotide identity.")

    doc.add_heading("S5. Complete MerMade 4 synthesis program", level=1)
    synthesis_rows = primary["tables"][4][1:]
    add_table(doc,
        ["Script detail", "Reagent / manufacturer", "No./V primes", "Injection", "Reaction", "Equalize"],
        synthesis_rows, [1.02, 2.62, 0.72, 0.72, 0.72, 0.8], font_size=6.65)

    doc.add_heading("S6. Primer record", level=1)
    primers = [
        ["SpAB F", "ATGGGTTCTTCTCACCACCACCA", "58.6", "123 A / 156 B"],
        ["SpA R", "TTAGCCGCAGTAGTTTTCCA", "52.3", "123 A"],
        ["SpB R", "TTAGCGGCGGGTTTTTGG", "53.8", "156 B"],
        ["SeqAB F", "TGCATCCATATGGGTTCTTCTCACCACCACCA", "67.2", "144 A / 177 B"],
        ["SeqA R", "TGCATCCTCGAGTTAGCCGCAGTAGTTTTCCA", "65.8", "144 A"],
        ["SeqB R", "TGCATCCTCGAGTTAGCGGCGGGTTTTTGG", "66.4", "177 B"],
        ["T7 F", "TAATACGACTCACTATAGGG", "52.0", "307 A / 340 B"],
        ["T7 R", "GCTAGTTATTGCTCAGCGG", "54.3", "307 A / 340 B"],
    ]
    add_table(doc, ["Primer", "Sequence (5′→3′)", "Reported Tm (°C)", "Expected product (bp)"], primers,
              [0.72, 3.4, 1.0, 1.48], font_size=7.8)
    add_text(doc,
        "The source table reported some A/B product sizes together in cells that do not apply to both reverse primers. The corrected product labels above associate SpA R/SeqA R with A and SpB R/SeqB R with B. Before a prospective run, each primer should be re-evaluated in G-Synth against the final recombinant plasmid and the laboratory's polymerase/buffer conditions."
    )

    doc.add_heading("S7. Wet-laboratory PCR and colony-screening evidence", level=1)
    add_figure(doc, SOURCE_FIG / "Figure_S3_Precloning_PCR.png",
        "Figure S3. PCR products obtained with insert-specific primers before cloning, reproduced from the primary experimental draft. SM, GeneRuler Express DNA Ladder; lanes 1–2, A insert; lanes 3–4, B insert.")
    add_figure(doc, SOURCE_FIG / "Figure_S4_Colony_PCR.png",
        "Figure S4. Colony-screening PCR reproduced from the primary experimental draft. Lanes 1 and 3 are A and B products with insert-specific primers; lanes 2 and 4 are the corresponding products with T7 primers.")

    doc.add_heading("S8. Sequencing evidence manifest", level=1)
    add_text(doc,
        "The validation set is limited to A Forward Seq.ab1, A Reverse Seq.ab1, B Forward Seq.ab1 and B Reverse Seq.ab1. Their SHA-256 hashes are recorded in the evidence JSON; no other chromatogram contributes to any reported result. Each F/R pair was oriented and assembled before consensus coverage was evaluated. Consensus coverage, identity, bidirectional overlap and F/R agreement are the primary retrospective endpoints.")
    rows = []
    for chain in ("A", "B"):
        d = seq["chains"][chain]
        consensus = d["gsynth_forward_reverse_consensus_raw_q0"]
        biopython_key = next(key for key in d if key.startswith("biopython_") and key.endswith("_local_alignment"))
        bio = d[biopython_key]
        rows.append([
            chain,
            str(d["reference"]["length_bp"]),
            f"{consensus['coverage_percent']:.1f}%",
            f"{consensus['identity_percent']:.1f}%",
            f"{consensus['bidirectional_overlap_percent']:.1f}%",
            f"{consensus['bidirectional_overlap_agreement_percent']:.1f}%",
            str(len(consensus["differences"])),
            f"{bio['combined_reference_coverage_percent']:.1f}%",
        ])
    add_table(doc,
        ["Chain", "Ref. bp", "Consensus cov.", "Identity", "F/R overlap", "F/R agree", "Consensus diff.", "Biopython cov."],
        rows, [0.45, 0.62, 0.86, 0.72, 0.78, 0.78, 0.82, 0.88], font_size=7.7)
    add_figure(doc, EVIDENCE / "sequencing_validation" / "Figure_GSynth_Sanger_Approved_Traces.png",
        "Figure S5. Primary G-Synth F/R consensus evidence map. The two oriented reads jointly cover each reference completely; green shading marks bidirectional overlap, with 100% agreement in both constructs.")

    doc.add_heading("S9. Supplementary interface evidence", level=1)
    add_figure(doc, INTERFACE / "Primer_Hybridization.png",
        "Figure S6. Primer–template hybridization view. The annealing 3′ segment is base-paired; the 5′ cloning extension is explicitly unpaired in cycle 1 and appears in the predicted product after extension.", width=4.0)
    add_figure(doc, INTERFACE / "Hybridization_Workflow.png",
        "Figure S7. Detailed SSD duplex hybridization. The two order molecules form a 123-bp exact antiparallel core; the NdeI-derived 5′-TA and XhoI-derived 5′-TCGA terminal products remain exposed and are labelled as cohesive ends.", width=6.35)
    add_figure(doc, INTERFACE / "Cloning_PreLigation.jpg",
        "Figure S8. Pre-ligation vector–insert compatibility gate. The simple view preserves enzyme identity, overhang sequence and physical polarity while withholding product analyses.", width=6.35)
    add_figure(doc, INTERFACE / "Cloning_Ligated_Product.jpg",
        "Figure S9. Product state after explicit in-silico ligation. The recombinant length and downstream map are displayed only after both junctions have passed.", width=6.35)
    add_figure(doc, INTERFACE / "HindIII_Diagnostic_Gel.png",
        "Figure S10. HindIII-containing diagnostic digest and predicted gel. The restriction site and complete digest fragments are displayed with a named broad-range marker. The panel is labelled as an in-silico prediction.", width=4.4)
    add_figure(doc, INTERFACE / "HindIII_Restriction_Map.png",
        "Figure S11. Restriction map showing HindIII as a named site derived from cut geometry.", width=5.4)

    doc.add_heading("S10. Secondary Geneious comparison and independent confirmation", level=1)
    add_text(doc,
        f"The supplied experimental draft used Benchling for reverse complements and Geneious Prime 2024.0.7 for cloning views and trace alignments. G-Synth independently regenerated the same four synthesis molecules from the two optimized coding inputs, and its recombinant maps placed the inserts downstream of the pET-21a(+) T7/lac expression elements with preserved coding frames. The archived Geneious figures and G-Synth therefore agree qualitatively on construct identity and placement. A new automated Geneious run was not performed because the reproducibility environment did not include a licensed command-line execution path. Biopython {biopython_version} was used as the independent, executable alignment comparator."
    )
    add_text(doc,
        "A prospective head-to-head study should freeze one set of reference plasmids and raw trace files, define identical trimming thresholds and variant rules, blind sample labels, and compare: read placement; orientation calls; covered intervals; mismatches/indels; consensus; runtime; manual interventions; exportability; and final verdict. The present study is a functional cross-check, not a statistical superiority benchmark."
    )
    add_figure(doc, SOURCE_FIG / "Figure_S6_Geneious_Duplex_Secondary_Comparison.png",
        "Figure S12. Archived Geneious Prime 2024.0.7 duplex representation of NdeI/XhoI-compatible A- and B-chain inserts. This is secondary qualitative concordance; G-Synth SSD is the primary in-silico reconstruction.")
    add_figure(doc, SOURCE_FIG / "Figure_S7_Geneious_Annotation_Secondary_Comparison.png",
        "Figure S13. Archived Geneious Prime construct annotation for the glargine A- and B-chain pET-21a(+) designs. G-Synth generated the primary coordinate, ORF and junction validation reported in the main article.")
    add_figure(doc, SOURCE_FIG / "Figure_S5_Geneious_Approved_Trace_Secondary_Comparison.png",
        "Figure S14. Archived Geneious Prime alignment made with the same four author-approved construct-confirmation files used by G-Synth. It provides secondary visual concordance and was not treated as an independent blinded benchmark.")

    doc.add_heading("S11. Scientific correction ledger for the supplied draft", level=1)
    add_table(doc,
        ["Source-draft implication", "Correction in the G-Synth manuscript", "Reason"],
        [
            ["Sequencing showed complete 100% identity", "Oriented F/R assembly produced 100% consensus coverage and identity; overlap support was 56.1% (A) and 71.2% (B), with 100% F/R agreement", "Separates assembled coverage from bidirectional confirmation without obscuring the complete consensus"],
            ["The method eliminates enzymatic digestion", "Insert preparation is PCR-free and pre-engineered; the pET-21a(+) vector is still digested with NdeI/XhoI", "Vector digestion is experimentally required"],
            ["The constructs enable insulin glargine production", "Construct design and cloning were validated; expression, cleavage, refolding and bioactivity were not tested", "DNA construction does not establish functional insulin"],
            ["His-tags and linkers ensure folding remains unaffected", "Sequence presence and tag position are reported; folding and cleavage remain experimental", "No protein-level evidence was supplied"],
            ["The workflow is broadly scalable and cost-effective", "Potential advantages are framed as hypotheses for prospective benchmarking", "No controlled cost, yield or scale comparison was performed"],
            ["Geneious is the primary validator", "G-Synth is the primary design/cloning/trace system; archived Geneious views and executable Biopython provide independent confirmation", "Matches the software paper's validation hierarchy"],
            ["One oligo length per insert", "Forward/reverse lengths are reported separately (125/127 nt and 158/160 nt)", "Cohesive-end geometry is strand-asymmetric"],
        ], [1.9, 2.9, 1.8], font_size=7.8)

    doc.add_heading("S12. Release and submission checklist", level=1)
    add_bullets(doc, [
        "Obtain written approval from every listed author for author order, affiliations, CRediT roles and AI disclosure.",
        "Run prospective bidirectional sequencing until both junctions and 100% of each insert are covered at the predefined quality threshold.",
        "Deposit raw trace files, reference FASTA/GenBank files, evidence JSON and analysis scripts in a stable repository, subject to institutional approval.",
        "Publish the version 1.0.0 GitHub release, connect Zenodo and insert the versioned DOI in the manuscript and CITATION.cff.",
        "Perform a clean-clone execution of all three test suites and both publication scripts; retain logs and environment versions.",
        "Prepare the target journal's graphical summary, high-resolution figures and machine-readable constructs; confirm its current author guidelines and retain the limitation that functional insulin glargine was not demonstrated.",
    ])

    PUBLICATION.mkdir(parents=True, exist_ok=True)
    doc.save(SI_OUT)


def main() -> None:
    design = json.loads(DESIGN_JSON.read_text())
    seq = json.loads(SEQ_JSON.read_text())
    primary = json.loads(PRIMARY_RECORD_JSON.read_text())
    build_main(design, seq, primary)
    build_si(design, seq, primary)
    print(MAIN_OUT)
    print(SI_OUT)


if __name__ == "__main__":
    main()
