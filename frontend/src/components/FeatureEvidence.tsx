export default function FeatureEvidence({ annotation }: {
  annotation: { name: string; type?: string; inferred?: boolean; basis?: string; regulatory_class?: string };
}) {
  const isRbs = annotation.type?.toLowerCase() === "rbs"
    || annotation.regulatory_class === "ribosome_binding_site" || /shine.dalgarno|sd-like/i.test(annotation.name);
  if (!annotation.inferred && !annotation.basis && !isRbs) return null;
  return <div className="feature-evidence">
    {annotation.inferred && <strong>Detected candidate · function unconfirmed</strong>}
    {annotation.basis && <p>{annotation.basis}</p>}
    {isRbs && <p>A Shine–Dalgarno motif can contribute to a bacterial RBS on mRNA, upstream of its translation initiation codon in the 5′→3′ direction. This map uses DNA coordinates. A motif match does not establish an active RBS or identify which protein it regulates; not all bacterial RBSs use a Shine–Dalgarno motif. <a href="https://pubmed.ncbi.nlm.nih.gov/38206950/" target="_blank" rel="noreferrer">Scientific reference</a></p>}
  </div>;
}
