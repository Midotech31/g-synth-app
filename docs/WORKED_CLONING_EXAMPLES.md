# Worked restriction-cloning examples

These examples are regression fixtures and teaching cases. Coordinates in the
interface are displayed one-based; engine slices remain zero-based internally.

## NdeI / XhoI: two 5′ cohesive ends

NdeI is a special coding-boundary case: cleavage retains an initiating ATG.
When “Use the site's ATG” is selected, PCR begins at the template's second
codon so the expressed product starts Met-X rather than Met-Met-X. XhoI leaves
a TCGA cohesive end at the 3′ boundary. A clean vector has one site for each;
recutting the recombinant plasmid must release exactly the inserted cassette.

## KpnI / SacI: two 3′ cohesive ends

The protruding strand is reversed relative to a 5′ cutter. Compatibility must
therefore be calculated from both strands and the side of the fragment, not
from a remembered overhang label. The clone/digest round trip is tested for
random inserts with this pair.

## NdeI / KpnI: mixed polarity

One junction is a 5′ cohesive end and the other is a 3′ cohesive end. The pair
forces orientation, but only if the insert ends are read from the actual duplex.
This case catches implementations that assume both ends share one polarity.

## EcoRV / SmaI: blunt / blunt

The molecules are compatible but orientation is not encoded by the ends and
background self-ligation is more likely. G-Synth treats identical enzyme use
as blocked and reports the experimental limitations of blunt-end ligation.

## Expected failure: internal recognition site

If the insert contains another selected recognition site, the PCR product or
recombinant plasmid cuts into additional fragments. Preflight reports
`PCR_RESTRICTION_SITE_COUNT` or the corresponding cloning gate and disables
release. Changing the enzyme pair or removing the internal site is required.
