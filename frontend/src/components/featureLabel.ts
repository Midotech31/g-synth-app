import type { Annotation } from "../api/client";


export function featureLabel(annotation: Pick<Annotation, "name" | "type" | "inferred">): string {
  if (!annotation.inferred) return annotation.name;
  if (/shine.dalgarno|sd-like/i.test(annotation.name)) {
    return annotation.type === "RBS" ? "? SD candidate" : "? SD motif";
  }
  return `? ${annotation.name}`;
}
