import { useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError } from "../api/client";
import { useAuth } from "../auth/AuthContext";

export default function Login() {
  const { signIn } = useAuth();
  const navigate = useNavigate();
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);

  async function onSubmit(event: React.FormEvent) {
    event.preventDefault();
    setError("");
    setBusy(true);
    try {
      await signIn(email, password);
      navigate("/", { replace: true });
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Could not sign in. Please try again.");
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="auth">
      <div className="panel">
        <div className="lockup">
          <div className="mark">🧬</div>
          <h1>G-Synth</h1>
          <p>Gene design workspace</p>
        </div>

        <div className="card">
          <div className="card-body">
            <form onSubmit={onSubmit}>
              <div className="field">
                <label htmlFor="email">Email</label>
                <input
                  id="email"
                  type="email"
                  autoComplete="email"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                  placeholder="you@lab.org"
                  required
                />
              </div>
              <div className="field">
                <label htmlFor="password">Password</label>
                <input
                  id="password"
                  type="password"
                  autoComplete="current-password"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                  required
                />
              </div>

              {error && <div className="notice notice-error">{error}</div>}

              <button className="btn btn-primary" type="submit" disabled={busy}>
                {busy && <span className="spinner" />}
                {busy ? "Signing in…" : "Sign in"}
              </button>
            </form>
          </div>
        </div>

        <p className="alt">
          No account yet? <Link to="/signup">Create one</Link>
        </p>
      </div>
    </div>
  );
}
