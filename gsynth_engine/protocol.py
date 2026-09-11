from __future__ import annotations

import csv
import io
from dataclasses import dataclass

from gsynth_engine.duplex import construct_duplex
from gsynth_engine.esd import ESDResult
from gsynth_engine.sequence import gc_content
from gsynth_engine.thermo import ANNEALING, melting_temperature


@dataclass(frozen=True)
class OligoOrder:


    name: str
    sequence: str
    length: int
    gc_percent: float
    tm: float
    scale: str
    purification: str
    role: str
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

    if length <= 60:
        return "25 nmol", "Desalted"
    if length <= 100:
        return "25 nmol", "PAGE"
    return "50 nmol", "PAGE"


def order_sheet(plan: ESDResult, *, construct_name: str = "construct") -> list[OligoOrder]:

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


def order_sheet_csv(plan: ESDResult, *, construct_name: str = "construct") -> str:

    orders = order_sheet(plan, construct_name=construct_name)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(orders[0].as_row) if orders else [])
    writer.writeheader()
    for order in orders:
        writer.writerow(order.as_row)  # type: ignore[arg-type]
    return buffer.getvalue()


def bench_protocol(
    plan: ESDResult,
    *,
    construct_name: str = "construct",
    vector: str = "pET-21a(+)",
) -> str:

    ssd = plan.ssd
    left, right = ssd.left_enzyme, ssd.right_enzyme
    n = plan.fragment_count
    orders = order_sheet(plan, construct_name=construct_name)

    lines: list[str] = []
    add = lines.append

    add(f"EXTENDED SEQUENCE DESIGN — {construct_name}")
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
        for left_frag, right_frag in zip(plan.fragments, plan.fragments[1:], strict=False):
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


def cloning_worksheet(
    result,
    *,
    vector_name: str,
    primer_set,
    ligation_plans: list,
    preflight,
    provenance: dict[str, object],
) -> str:

    lines: list[str] = []
    add = lines.append
    add(f"G-SYNTH CLONING WORKSHEET — {result.name}")
    add("=" * 76)
    add(f"Vector             {vector_name}")
    add(f"Enzyme pair        {result.left_enzyme} / {result.right_enzyme}")
    add(f"Recombinant        {result.length} bp")
    add(f"Backbone band      {result.backbone_length} bp after double digest")
    add(f"Excised vector     {result.removed_length} bp")
    add(f"Insert             {result.insert_length} bp")
    add(f"Engine             {provenance.get('engine_version', 'unknown')}")
    add(f"Output SHA-256     {provenance.get('output_sha256', 'unavailable')}")
    add(f"Parameters SHA-256 {provenance.get('parameters_sha256', 'unavailable')}")
    add("")

    add("1. PREFLIGHT RELEASE")
    add("-" * 76)
    add(f"Overall verdict    {preflight.verdict.upper()}")
    for check in preflight.checks:
        mark = "PASS" if check.status == "pass" else check.status.upper()
        add(f"[{mark:6}] {check.code:<28} {check.label}")
        add(f"         {check.detail}")
    add("")

    add("2. DOUBLE-DIGEST RECORD")
    add("-" * 76)
    add("Confirm buffer, enzyme units, temperature and heat inactivation against")
    add("the current supplier datasheets; record the lot-specific setup below.")
    add(f"Expected gel bands: {result.backbone_length} bp backbone + "
        f"{result.removed_length} bp excised vector segment.")
    add("Buffer/lot: ____________________  Incubation: ____________________")
    add("Observed bands: __________________  Purified backbone: __________ ng/µL")
    add("")

    add("3. LIGATION SERIES")
    add("-" * 76)
    add("Use the measured purified-DNA concentrations to convert these target")
    add("masses to pipetting volumes. Include vector-only and no-ligase controls.")
    for plan in ligation_plans:
        add(f"{plan.ratio:g}:1 insert:vector — vector {plan.vector_ng:g} ng "
            f"({plan.vector_fmol:g} fmol), insert {plan.insert_ng:g} ng "
            f"({plan.insert_fmol:g} fmol)")
    add("Vector-only control: ______  No-ligase control: ______  Plate ID: ______")
    add("")

    add("4. SEQUENCING PRIMERS")
    add("-" * 76)
    for primer in primer_set.primers:
        direction = "forward" if primer.direction == 1 else "reverse"
        add(f"{primer.name:<24} {primer.sequence}  {primer.tm:.1f}°C  {direction}")
        add(f"{'':24} expected read {primer.reads_from + 1}–{primer.reads_to}")
    if primer_set.gaps:
        add(f"COVERAGE WARNING: uncovered intervals {primer_set.gaps}")
    else:
        add("Coverage gate: the complete insert is covered by the primer set.")
    add("")

    add("5. RESULTS AND SIGN-OFF")
    add("-" * 76)
    add("Colony/clone ID: __________________  Extraction date: ________________")
    add("Digest result:  □ expected  □ unexpected  □ not run")
    add("Sequencing:     □ fully verified  □ differences  □ partial/unplaced")
    add("Reviewed by: ______________________  Date: __________________________")
    add("")
    add("Do not release a clone for expression until the requested region is")
    add("fully covered and every confident difference has been adjudicated.")
    return "\n".join(lines) + "\n"
