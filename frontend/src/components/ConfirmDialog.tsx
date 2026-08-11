import { useCallback, useEffect, useId, useRef } from "react";

/**
 * The question asked before something irreversible happens.
 *
 * This replaced `window.confirm`, which is announced by the browser rather
 * than by the page: it cannot say *which* project is about to go, it stops
 * every other script while it waits, and on some screen readers it is read
 * out with no indication that the surrounding page is now inert. Deleting a
 * design is not recoverable from the interface, so the one moment that
 * stands between a person and losing it should be the part that behaves.
 *
 * Focus moves in, is held inside while the question is open, and goes back
 * to the button that asked it — a keyboard user who declines ends up where
 * they were rather than at the top of the document.
 */

type Props = {
  open: boolean;
  /** Names the thing, not the operation: "Delete “pGS-EntA”?" */
  title: string;
  /** What happens if they say yes, and whether it can be undone. */
  body: React.ReactNode;
  confirmLabel?: string;
  cancelLabel?: string;
  onConfirm: () => void;
  onCancel: () => void;
};

const FOCUSABLE =
  'a[href], button:not([disabled]), textarea:not([disabled]), input:not([disabled]), ' +
  'select:not([disabled]), [tabindex]:not([tabindex="-1"])';

export default function ConfirmDialog({
  open,
  title,
  body,
  confirmLabel = "Delete",
  cancelLabel = "Cancel",
  onConfirm,
  onCancel,
}: Props) {
  const panel = useRef<HTMLDivElement>(null);
  const declineButton = useRef<HTMLButtonElement>(null);
  const opener = useRef<HTMLElement | null>(null);
  const titleId = useId();
  const bodyId = useId();

  useEffect(() => {
    if (!open) return;

    // Whatever had focus is where the question was asked from, and where the
    // answer has to leave the user. A trigger that the answer itself removes
    // from the page — the delete button on the card being deleted — cannot
    // be returned to; the caller moves focus somewhere that still exists.
    opener.current = document.activeElement as HTMLElement | null;
    // The declining option takes focus, not the destructive one: the cost of
    // a stray Enter here is a design that no longer exists.
    declineButton.current?.focus();

    return () => {
      const previous = opener.current;
      if (previous && previous.isConnected) previous.focus();
    };
  }, [open]);

  const onKeyDown = useCallback(
    (event: React.KeyboardEvent) => {
      if (event.key === "Escape") {
        event.stopPropagation();
        onCancel();
        return;
      }
      if (event.key !== "Tab") return;

      // Tab off either end wraps, so the rest of the page cannot be reached
      // while a question about it is unanswered.
      const stops = Array.from(
        panel.current?.querySelectorAll<HTMLElement>(FOCUSABLE) ?? [],
      );
      if (stops.length === 0) return;
      const first = stops[0];
      const last = stops[stops.length - 1];
      const active = document.activeElement;

      if (event.shiftKey && (active === first || active === panel.current)) {
        event.preventDefault();
        last.focus();
      } else if (!event.shiftKey && active === last) {
        event.preventDefault();
        first.focus();
      }
    },
    [onCancel],
  );

  if (!open) return null;

  return (
    <div
      className="modal-backdrop"
      // Clicking away is the same answer as Cancel, but only when the click
      // began outside: a drag that starts on the text must not dismiss it.
      onMouseDown={(event) => {
        if (event.target === event.currentTarget) onCancel();
      }}
    >
      <div
        className="modal"
        ref={panel}
        role="dialog"
        aria-modal="true"
        aria-labelledby={titleId}
        aria-describedby={bodyId}
        tabIndex={-1}
        onKeyDown={onKeyDown}
      >
        <h2 id={titleId}>{title}</h2>
        <p id={bodyId}>{body}</p>
        <div className="modal-actions">
          <button type="button" className="btn btn-outline" ref={declineButton} onClick={onCancel}>
            {cancelLabel}
          </button>
          <button type="button" className="btn btn-danger-solid" onClick={onConfirm}>
            {confirmLabel}
          </button>
        </div>
      </div>
    </div>
  );
}
