/** Display evidence and biological terminology consistently across viewers. */
export default function FeatureEvidence({ annotation }: {
  annotation: { name: string; type?: string; inferred?: boolean; basis?: string; regulatory_class?: string };
}) {
  const isRbs = annotation.type?.toLowerCase() === "rbs"
    || annotation.regulatory_class === "ribosome_binding_site" || /shine.dalgarno/i.test(annotation.name);
  if (!annotation.inferred && !annotation.basis && !isRbs) return null;
  return <div className="feature-evidence">
    {annotation.inferred && <strong>Detected candidate · function unconfirmed</strong>}
    {annotation.basis && <p>{annotation.basis}</p>}
    {isRbs && <p>The Shine–Dalgarno motif acts on bacterial mRNA, upstream of the translation initiation codon in the 5′→3′ reading direction. This DNA map marks the corresponding transcribed region; a motif match alone does not establish an active RBS. <a href="https://pubmed.ncbi.nlm.nih.gov/38206950/" target="_blank" rel="noreferrer">Scientific reference</a></p>}
  </div>;
}
