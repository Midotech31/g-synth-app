"""From a design to the bench: an order sheet and a working protocol.

This is where G-Synth goes further than a sequence viewer. The design is not
the deliverable — the deliverable is the list of oligos to order and the
steps to run once they arrive.

Nothing here invents numbers. Volumes and temperatures are the standard
conditions for annealing synthetic oligos and ligating them with T4 DNA
ligase; anything design-specific (oligo count, overhangs, fragment order)
comes from the :class:`AssemblyPlan`.
"""
from __future__ import annotations

import csv
import io
from dataclasses import dataclass

from gsynth_engine.duplex import construct_duplex
from gsynth_engine.merzoug import AssemblyPlan
from gsynth_engine.sequence import gc_content
from gsynth_engine.thermo import ANNEALING, melting_temperature


@dataclass(frozen=True)
class OligoOrder:
    """One line of the order sheet."""

    name: str
    sequence: str
    length: int
    gc_percent: float
    tm: float
    scale: str
    purification: str
    role: str          # "forward" | "reverse"
    fragment: int

    @property
    def as_row(self) -> dict[str, object]:
        return {
            "Name": self.name,
            "Sequence (5'->3')": self.sequence,
            "Length (nt)": self.length,
            "GC (%)": self.gc_percent,
            "Tm (°C)": self.tm,
            "Scale": self.scale,
            "Purification": self.purification,
            "Fragment": self.fragment,
            "Strand": self.role,
        }


def _recommend_scale(length: int) -> tuple[str, str]:
    """Synthesis scale and purification most suppliers would advise.

    Longer oligos accumulate more truncated products, so past ~60 nt a
    purification step stops those failures ending up in the ligation.
    """
    if length <= 60:
        return "25 nmol", "Desalted"
    if length <= 100:
        return "25 nmol", "PAGE"
    return "50 nmol", "PAGE"


def order_sheet(plan: AssemblyPlan, *, construct_name: str = "construct") -> list[OligoOrder]:
    """Every oligo to order, in the order fragments will be assembled."""
    prefix = construct_name.strip().replace(" ", "_") or "construct"
    orders: list[OligoOrder] = []
    for fragment in plan.fragments:
        for role, sequence in (("forward", fragment.forward), ("reverse", fragment.reverse)):
            scale, purification = _recommend_scale(len(sequence))
            suffix = "F" if role == "forward" else "R"
            orders.append(
                OligoOrder(
                    name=f"{prefix}_{fragment.name}_{suffix}",
                    sequence=sequence,
                    length=len(sequence),
                    gc_percent=round(gc_content(sequence), 1),
                    tm=round(melting_temperature(sequence, conditions=ANNEALING), 1),
                    scale=scale,
                    purification=purification,
                    role=role,
                    fragment=fragment.index,
                )
            )
    return orders


def order_sheet_csv(plan: AssemblyPlan, *, construct_name: str = "construct") -> str:
    """The order sheet as CSV — most suppliers accept an upload in this shape."""
    orders = order_sheet(plan, construct_name=construct_name)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(orders[0].as_row) if orders else [])
    writer.writeheader()
    for order in orders:
        writer.writerow(order.as_row)  # type: ignore[arg-type]
    return buffer.getvalue()


def bench_protocol(
    plan: AssemblyPlan,
    *,
    construct_name: str = "construct",
    vector: str = "pET-21a(+)",
) -> str:
    """A protocol a student can follow without reading the source code.

    Covers resuspension, annealing, the pairwise ligation Merzoug assembly
    calls for, and cloning into the cut vector.
    """
    ssd = plan.ssd
    left, right = ssd.left_enzyme, ssd.right_enzyme
    n = plan.fragment_count
    orders = order_sheet(plan, construct_name=construct_name)

    lines: list[str] = []
    add = lines.append

    add(f"MERZOUG ASSEMBLY — {construct_name}")
    add("=" * 72)
    add("")
    add(f"Construct        {plan.construct_length} bp "
        f"({gc_content(plan.construct_forward):.1f}% GC)")
    add(f"Fragments        {n}  →  {plan.oligo_count} oligos to order")
    add(f"Junctions        {len(plan.junction_overhangs)} × "
        f"{plan.overhang_length} nt 5' overhangs: "
        + (", ".join(plan.junction_overhangs) or "none"))
    add(f"Cloning ends     {left} (5'-{ssd.left_overhang}) / "
        f"{right} (5'-{ssd.right_overhang})")
    add(f"Vector           {vector}, cut with {left} + {right}")
    if ssd.cleavage_site:
        add(f"Tag              6×His + {ssd.cleavage_site} site")
    add("")
    add("No PCR is used at any step. Internal junctions are ligated through")
    add("synthetic overhangs — no restriction digestion of the fragments.")
    add("")

    add("1. OLIGOS TO ORDER")
    add("-" * 72)
    width = max(len(o.name) for o in orders)
    for order in orders:
        add(f"   {order.name:<{width}}  {order.length:>3} nt  "
            f"Tm {order.tm:>5.1f}°C  {order.scale:>9}  {order.purification}")
    add("")
    add("   Tm: nearest-neighbour model (SantaLucia 1998) under the annealing")
    add(f"   conditions of step 3 — {ANNEALING.summary}.")
    add("")

    add("2. RESUSPENSION")
    add("-" * 72)
    add("   Resuspend each oligo to 100 µM in TE (pH 8.0). Vortex, spin down,")
    add("   and leave 15 min at room temperature before use.")
    add("")

    add("3. ANNEALING — one reaction per fragment")
    add("-" * 72)
    add("   Per fragment, mix:")
    add("       forward oligo (100 µM)        5 µL")
    add("       reverse oligo (100 µM)        5 µL")
    add("       10× annealing buffer          2 µL")
    add("       nuclease-free water           8 µL")
    add("                                    ------")
    add("                                     20 µL   (25 µM duplex)")
    add("")
    add("   Heat to 95 °C for 5 min in a heat block, then switch the block")
    add("   off and let it cool to room temperature over ~1 h. Slow cooling")
    add("   matters: a fast drop traps mispaired oligos.")
    add("")
    for fragment in plan.fragments:
        ends = (f"left 5'-{fragment.left_overhang}"
                f"{' (vector)' if fragment.is_first else ''}"
                f" · right 5'-{fragment.right_overhang}"
                f"{' (vector)' if fragment.is_last else ''}")
        add(f"       {fragment.name}: {len(fragment.forward)} + "
            f"{len(fragment.reverse)} nt — {ends}")
    add("")

    # The molecule itself, not a summary of it. Check this before ordering:
    # it is the only place a wrong overhang shows as a wrong overhang.
    add("   HYBRIDISATION — check this before the oligos are ordered")
    add("   " + "-" * 68)
    view = construct_duplex(plan)
    if view.mismatches():
        add(f"   !! {len(view.mismatches())} positions do not pair. Do not order.")
        add("")
    for line in view.to_text(60).rstrip().splitlines():
        add(f"   {line}" if line else "")
    add("")
    add("   Unpaired bases are single-stranded: the sticky ends at the two")
    add("   outer ends, and the junction overhangs where each strand is cut")
    add(f"   at a different position ({plan.overhang_length} nt apart).")
    add("")

    add("4. PHOSPHORYLATION")
    add("-" * 72)
    add("   Synthetic oligos carry no 5' phosphate, so ligase cannot join")
    add("   them. Phosphorylate every annealed duplex EXCEPT the two outer")
    add("   ends that go into the vector — leaving those unphosphorylated")
    add("   suppresses vector self-ligation.")
    add("       annealed duplex              10 µL")
    add("       10× T4 PNK buffer             2 µL")
    add("       10 mM ATP                     2 µL")
    add("       T4 polynucleotide kinase      1 µL")
    add("       water                         5 µL")
    add("   37 °C, 30 min, then 65 °C, 20 min to inactivate.")
    add("")

    add("5. LIGATION — pairwise, in order")
    add("-" * 72)
    if n == 1:
        add("   Single fragment: go straight to step 6.")
    else:
        add("   Combine fragments in equimolar amounts. Each junction has its")
        add("   own overhang, so the order below is enforced by the sequence")
        add("   itself — but ligating pairwise keeps the yield of correct")
        add("   product highest:")
        add("")
        chain = " + ".join(f.name for f in plan.fragments)
        add(f"       {chain}")
        add("")
        for left_frag, right_frag in zip(plan.fragments, plan.fragments[1:]):
            add(f"       {left_frag.name} → {right_frag.name} "
                f"via 5'-{left_frag.right_overhang}")
        add("")
        add("       each duplex (25 µM)         1 µL")
        add("       10× T4 ligase buffer        2 µL")
        add("       T4 DNA ligase               1 µL")
        add("       water                       to 20 µL")
        add("   16 °C overnight, or room temperature for 2 h.")
    add("")

    add("6. CLONING")
    add("-" * 72)
    add(f"   Digest {vector} with {left} + {right}, dephosphorylate the")
    add("   backbone (CIP, 37 °C, 30 min), and gel-purify it.")
    add("")
    add("   Ligate the assembled insert into the cut vector at a 3:1")
    add("   insert:vector molar ratio, 16 °C overnight.")
    add("")
    add("   Transform into a cloning strain (DH5α), select on the vector's")
    add("   antibiotic, and screen colonies by colony PCR or restriction")
    add("   digest. Sequence-verify before expressing.")
    add("")

    if ssd.warnings or plan.warnings:
        add("NOTES AND WARNINGS")
        add("-" * 72)
        for warning in dict.fromkeys([*ssd.warnings, *plan.warnings]):
            add(f"   • {warning}")
        add("")

    add("7. EXPECTED CONSTRUCT")
    add("-" * 72)
    add("   Forward strand (5'->3'):")
    for i in range(0, len(plan.construct_forward), 60):
        add(f"   {i + 1:>6}  {plan.construct_forward[i : i + 60]}")
    add("")
    add("   Verify by sequencing across every junction before expressing.")

    return "\n".join(lines)
