import { Component, type ErrorInfo, type ReactNode } from "react";

/**
 * The last thing between a rendering fault and a blank page.
 *
 * React unmounts the whole tree when a render throws, so without this one bad
 * value — a NaN coordinate reaching the plasmid map, an undefined array from a
 * response shape that changed — takes the entire workspace with it, including
 * whatever design was on screen and unsaved. What the reader sees instead is a
 * page that says which part failed and offers the two things that actually
 * recover: try the view again, or go back to the workspace.
 *
 * Deliberately a class. Error boundaries have no hook equivalent — there is no
 * `useErrorBoundary`, and `componentDidCatch` is the only way to intercept a
 * descendant's render error.
 */

type Props = {
  children: ReactNode;
  /** Names the part that failed, so the message is specific. */
  area?: string;
};

type State = { error: Error | null };

export default class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo): void {
    // Kept: the stack is the only record of what happened, and a lab machine
    // is not going to have a session replay to consult afterwards.
    console.error("Render failed:", error, info.componentStack);
  }

  private reset = (): void => {
    this.setState({ error: null });
  };

  render(): ReactNode {
    const { error } = this.state;
    if (!error) return this.props.children;

    const area = this.props.area ?? "this page";

    return (
      <div className="content">
        <div className="card" role="alert">
          <div className="card-head">
            <h2>Something went wrong in {area}</h2>
          </div>
          <div className="card-body" style={{ display: "flex", flexDirection: "column", gap: "0.9rem" }}>
            <p className="note">
              The page stopped before it finished drawing, so what you can see may be
              incomplete. Nothing you have saved is affected &mdash; results already
              downloaded or stored in a project are untouched.
            </p>
            <p className="note">
              Try the view again. If it fails the same way a second time, the input that
              caused it is likely still in the form, so change it before retrying.
            </p>

            <details>
              <summary className="note" style={{ cursor: "pointer" }}>
                Technical detail
              </summary>
              <pre className="seq-block" style={{ marginTop: "0.5rem", whiteSpace: "pre-wrap" }}>
                {error.message || String(error)}
              </pre>
            </details>

            <div style={{ display: "flex", gap: "0.5rem", flexWrap: "wrap" }}>
              <button className="btn btn-primary" onClick={this.reset}>
                Try again
              </button>
              {/* A full load, not a route change: the fault may have left state
                  that a re-render would walk straight back into. */}
              <button className="btn btn-outline" onClick={() => { window.location.href = "/"; }}>
                Back to workspace
              </button>
            </div>
          </div>
        </div>
      </div>
    );
  }
}
