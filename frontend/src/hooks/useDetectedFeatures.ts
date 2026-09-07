import { useEffect, useMemo, useState } from "react";
import { api, type Annotation, type DetectedFeature } from "../api/client";

/** Preview candidates automatically; only explicit acceptance changes saved features. */
export function useDetectedFeatures(sequence: string, annotations: Annotation[], circular: boolean) {
  const [scan, setScan] = useState<{ sequence: string; annotations: Annotation[]; matches: DetectedFeature[]; error: string } | null>(null);
  const [revision, setRevision] = useState(0);
  useEffect(() => {
    if (!sequence) return;
    let cancelled = false;
    const timer = setTimeout(() => {
      api.detectSequenceFeatures(sequence, annotations, circular).then((result) => {
        if (!cancelled) setScan({ sequence, annotations, matches: result.matches, error: "" });
      }).catch(() => {
        if (!cancelled) setScan({ sequence, annotations, matches: [], error: "Automatic feature detection is unavailable. Retry the scan." });
      });
    }, 200);
    return () => { cancelled = true; clearTimeout(timer); };
  }, [sequence, annotations, circular, revision]);
  const current = scan?.sequence === sequence && scan.annotations === annotations ? scan : null;
  const matches = useMemo(() => current?.matches ?? [], [current]);
  return { matches, scanning: !!sequence && !current, error: current?.error ?? "", rescan: () => { setScan(null); setRevision((value) => value + 1); } };
}
