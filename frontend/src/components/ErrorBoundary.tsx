import { Component, type ErrorInfo, type ReactNode } from "react";


type Props = {
  children: ReactNode;

  area?: string;
};

type State = { error: Error | null };

export default class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo): void {


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
