import { useEffect, useRef, useState, type ReactNode } from "react";

const FOCUSABLE = 'button:not([disabled]), a[href], input:not([disabled]), select:not([disabled]), textarea:not([disabled]), [tabindex="0"]';

/** Expand the existing viewer without losing its sequence, zoom or selection. */
export default function ExpandablePanel({ children, label, className = "" }: {
  children: ReactNode; label: string; className?: string;
}) {
  const [expanded, setExpanded] = useState(false);
  const panel = useRef<HTMLElement>(null);
  const toggle = useRef<HTMLButtonElement>(null);
  useEffect(() => {
    if (!expanded) return;
    const previous = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    toggle.current?.focus();
    return () => { document.body.style.overflow = previous; toggle.current?.focus(); };
  }, [expanded]);
  return (
    <section ref={panel} className={`expandable-panel ${className}${expanded ? " expanded" : ""}`}
      role={expanded ? "dialog" : "region"} aria-modal={expanded || undefined} aria-label={label}
      onKeyDown={(event) => {
        if (!expanded) return;
        if (event.key === "Escape") { event.stopPropagation(); setExpanded(false); }
        if (event.key !== "Tab") return;
        const stops = [...(panel.current?.querySelectorAll<HTMLElement>(FOCUSABLE) ?? [])]
          .filter((element) => !element.closest('[hidden]'));
        const first = stops[0]; const last = stops[stops.length - 1];
        if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last?.focus(); }
        if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first?.focus(); }
      }}>
      <div className="viewer-expand-toolbar">
        <span>{label}</span>
        <button ref={toggle} type="button" className="btn btn-outline" aria-expanded={expanded}
          onClick={() => setExpanded(!expanded)}>{expanded ? "Restore viewer" : "Expand viewer"}</button>
      </div>
      {children}
    </section>
  );
}
