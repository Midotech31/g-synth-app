import { useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError } from "../api/client";
import { useAuth } from "../auth/AuthContext";
import AuthShell from "../components/AuthShell";
import PasswordInput from "../components/PasswordInput";

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
    <AuthShell
      eyebrow="Secure workspace"
      title="Welcome back"
      description="Continue your design, build and validation work from one scientific record."
      footer={<>No account yet? <Link to="/signup">Create one</Link></>}
    >
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
                <PasswordInput
                  id="password"
                  label="Password"
                  autoComplete="current-password"
                  value={password}
                  onChange={setPassword}
                />
              </div>

              {/* A failed sign-in is the one thing on this page a reader must
                  be told about, and it appears without the focus moving. */}
              {error && <div className="notice notice-error" role="alert">{error}</div>}

              <button className="btn btn-primary" type="submit" disabled={busy}>
                {busy && <span className="spinner" />}
                {busy ? "Signing in…" : "Sign in"}
              </button>
      </form>
    </AuthShell>
  );
}
