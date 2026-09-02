#!/usr/bin/env python3
"""Build deterministic publication figures from G-Synth evidence and UI captures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "publication_evidence" / "manuscript_figures"
INTERFACE = ROOT / "publication_evidence" / "interface_evidence"


def _workflow() -> None:
    fig, ax = plt.subplots(figsize=(13, 7.2), dpi=300)
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 7.2)
    ax.axis("off")
    navy, teal, green, purple, orange = "#12233F", "#39757D", "#4C8B66", "#705493", "#B06D36"
    colors = [teal, purple, orange, green]
    titles = ["DESIGN", "BUILD", "CLONE", "VERIFY"]
    subtitles = [
        "Codon optimisation\nrestriction-aware cassette\norderable oligos",
        "SSD / ESD reconstruction\norderable oligo pairs\nPCR / gel simulation",
        "Vector digestion\ndirectional ligation\nannotation + ORF checks",
        "AB1/SCF traces\nquality-gated alignment\nmismatch reporting",
    ]
    x_positions = [0.6, 3.75, 6.9, 10.05]
    for i, (x, c, title, subtitle) in enumerate(zip(x_positions, colors, titles, subtitles, strict=True)):
        box = FancyBboxPatch((x, 3.65), 2.35, 2.15, boxstyle="round,pad=0.12,rounding_size=0.18",
                             linewidth=1.5, edgecolor=c, facecolor="#F8FAFC")
        ax.add_patch(box)
        ax.text(x + 1.175, 5.25, title, ha="center", va="center", fontsize=14, fontweight="bold", color=c)
        ax.text(x + 1.175, 4.45, subtitle, ha="center", va="center", fontsize=10.5, color=navy, linespacing=1.35)
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + 2.42, 4.73), (x_positions[i + 1] - 0.08, 4.73),
                                         arrowstyle="-|>", mutation_scale=16, linewidth=2, color="#78909C"))

    ax.text(6.5, 6.65, "G-Synth: one auditable design-to-sequence record", ha="center",
            fontsize=21, fontweight="bold", color=navy)
    ax.text(6.5, 6.25, "Deterministic scientific calculations; editable biological context; explicit confidence boundaries",
            ha="center", fontsize=11.5, color="#506477")

    layers = [
        (0.8, 2.35, 11.4, 0.64, "React + TypeScript workspace", "#E8F1F3", teal),
        (0.8, 1.52, 11.4, 0.64, "Django REST API, accounts and project provenance", "#EEF0F8", purple),
        (0.8, 0.69, 11.4, 0.64, "Dependency-free Python scientific engine + 1,113 tests", "#EAF3EC", green),
    ]
    for x, y, w, h, label, fill, edge in layers:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.12",
                                    linewidth=1.1, edgecolor=edge, facecolor=fill))
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=11.2,
                fontweight="bold", color=navy)
    ax.text(0.22, 1.85, "SOFTWARE\nLAYERS", ha="center", va="center", rotation=90,
            fontsize=9.5, fontweight="bold", color="#657786")
    fig.tight_layout(pad=0.25)
    fig.savefig(OUT / "Figure_1_GSynth_workflow_and_architecture.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _font(size: int, bold: bool = False):
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ]
    for p in candidates:
        if Path(p).exists():
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def _fit(img: Image.Image, width: int) -> Image.Image:
    ratio = width / img.width
    return img.resize((width, round(img.height * ratio)), Image.Resampling.LANCZOS)


def _panel(paths: list[Path], labels: list[str], out_name: str, width: int = 1550, gap: int = 36) -> None:
    images = [_fit(Image.open(p).convert("RGB"), width) for p in paths]
    label_h = 74
    total_h = sum(im.height + label_h for im in images) + gap * (len(images) - 1)
    canvas = Image.new("RGB", (width, total_h), "white")
    draw = ImageDraw.Draw(canvas)
    y = 0
    for im, label in zip(images, labels, strict=True):
        draw.text((8, y + 10), label, fill="#12233F", font=_font(36, bold=True))
        y += label_h
        canvas.paste(im, (0, y))
        y += im.height + gap
    canvas.save(OUT / out_name, quality=95)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _workflow()
    _panel(
        [INTERFACE / "Annotated_Glargine_A.jpg", INTERFACE / "Annotated_Glargine_B.jpg"],
        ["A  Insulin glargine A-chain construct", "B  Insulin glargine B-chain construct"],
        "Figure_3_GSynth_annotated_glargine_constructs.jpg",
    )
    _panel(
        [INTERFACE / "Reference_Alignment_A.png", INTERFACE / "Reference_Alignment_B.png"],
        ["A  Reference-aligned A-chain chromatograms", "B  Reference-aligned B-chain chromatograms"],
        "Figure_4_GSynth_reference_aligned_viewer.jpg",
        width=1800,
    )
    print(f"Wrote manuscript figures to {OUT}")


if __name__ == "__main__":
    main()
