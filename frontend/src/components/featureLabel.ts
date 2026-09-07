import type { Annotation } from "../api/client";

/** Put uncertainty first so it survives clipping on short sequence features. */
export function featureLabel(annotation: Pick<Annotation, "name" | "type" | "inferred">): string {
  if (!annotation.inferred) return annotation.name;
  if (/shine.dalgarno|sd-like/i.test(annotation.name)) {
    return annotation.type === "RBS" ? "? SD candidate" : "? SD motif";
  }
  return `? ${annotation.name}`;
}
