# Vector reference checks

The user-provided MilliporeSigma (Novagen) references were checked on
2026-09-07 using the downloads linked from these SnapGene pages:

| Vector | Reference | Length | Comparison result |
| --- | --- | ---: | --- |
| pET-21a(+) | [SnapGene](https://www.snapgene.com/plasmids/pet_and_duet_vectors_(novagen)/pET-21a(%2B)) | 5,443 bp | Exact identity with the existing laboratory reference, including origin and orientation. Existing sequence and annotations are unchanged. |
| pET-28a(+) | [SnapGene](https://www.snapgene.com/plasmids/pet_and_duet_vectors_(novagen)/pET-28a(%2B)) | 5,369 bp | Matches the catalogue's expected length and unique restriction sites. Import remains required. |

The uppercase DNA SHA-256 fingerprints are:

- pET-21a(+): `7319006d8801ea36b8fdc5dbe508ed7f5ca4f0ccae2a8eb96966f4c2205d7ff0`
- pET-28a(+): `00f6f8d5571b3f503417f87feb2cd8ed47700f3201be18e009366ea793f8901d`

No new SnapGene sequence files, maps, annotations or descriptive notes are
published with this reference update. Source attribution:
[SnapGene resources](https://www.snapgene.com/resources).

For the complete G-Synth insert cassette, the pET-28a catalogue presents
NcoI / XhoI first because this replaces the native N-terminal leader. NdeI
and BamHI retain upstream leader sequence. The current insert-origin
translation does not certify that complete native fusion; review of the
expression context remains necessary for that case.
