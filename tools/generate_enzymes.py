from __future__ import annotations

import sys
from collections import defaultdict

try:
    from Bio import Restriction
except ImportError:
    sys.exit("Biopython is needed to regenerate this table: pip install biopython")

sys.path.insert(0, ".")
from gsynth_engine.constants import RESTRICTION_ENZYMES as CURATED  # noqa: E402


def usable(enzyme) -> bool:

    site = str(enzyme.site)
    if not site or set(site) - set("ACGT"):
        return False
    if enzyme.fst5 is None or enzyme.fst3 is None:
        return False
    if enzyme.cut_twice():
        return False
    top, bottom = enzyme.fst5, len(site) + enzyme.fst3
    return 0 <= top <= len(site) and 0 <= bottom <= len(site)


def canonical(names: list[str]) -> str:

    for name in names:
        if name in CURATED:
            return name
    return min(names, key=lambda n: (-len(getattr(Restriction, n).supplier_list()),
                                     len(n), n))


groups: dict[tuple[str, int, int], list[str]] = defaultdict(list)
for enzyme in Restriction.CommOnly:
    if usable(enzyme):
        site = str(enzyme.site)
        groups[(site, enzyme.fst5, len(site) + enzyme.fst3)].append(str(enzyme))

rows = []
for (site, top, bottom), names in groups.items():
    name = canonical(sorted(names))
    aliases = tuple(sorted(n for n in names if n != name))
    rows.append((name, site, top, bottom, aliases))
rows.sort()

print('"""Restriction-enzyme cut geometries used for detection and design.')
print()
print("Versioned REBASE-derived cut geometries processed with Biopython.")
print("Each entry is a distinct cut specification; isoschizomers are retained")
print("as aliases.")
print()
print("The preferred subset in ``constants.py`` is offered first in selectors.")
print('"""')
print("from __future__ import annotations")
print()
print("from typing import Final")
print()
print(f"#: {len(rows)} distinct specifications, from "
      f"{sum(len(v) for v in groups.values())} commercially available enzymes.")
print("ENZYME_TABLE: Final[dict[str, dict[str, object]]] = {")
for name, site, top, bottom, aliases in rows:
    alias = f', "aliases": {aliases!r}' if aliases else ""
    print(f'    "{name}": {{"recognition": "{site}", '
          f'"cut_top": {top}, "cut_bottom": {bottom}{alias}}},')
print("}")
