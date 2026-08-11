/**
 * What just happened, in one sentence, for someone who cannot see it happen.
 *
 * Every workspace here is the same shape: press the button, wait, read the
 * answer. Sighted users watch the spinner stop and the cards appear; nothing
 * about that moves focus, so a screen reader says nothing at all — the page
 * simply goes quiet and stays quiet. This says the verdict instead, because
 * "verified" and "do not order" are the whole point of the wait.
 *
 * Render it unconditionally: a live region that arrives on the page at the
 * same moment as its text is not announced, only one already there whose
 * contents then change.
 */
export default function LiveStatus({ message }: { message: string }) {
  return (
    <p className="sr-only" role="status" aria-live="polite">
      {message}
    </p>
  );
}
