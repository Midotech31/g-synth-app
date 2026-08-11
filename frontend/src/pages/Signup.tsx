import { useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError } from "../api/client";
import { useAuth } from "../auth/AuthContext";
import { Logo } from "../components/Logo";

export default function Signup() {
  const { signUp } = useAuth();
  const navigate = useNavigate();
  const [name, setName] = useState("");
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [password2, setPassword2] = useState("");
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);

  async function onSubmit(event: React.FormEvent) {
    event.preventDefault();
    setError("");
    setBusy(true);
    try {
      await signUp(email, name, password, password2);
      navigate("/", { replace: true });
    } catch (err) {
      setError(
        err instanceof ApiError ? err.message : "Could not create the account. Please try again.",
      );
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="auth">
      <div className="panel">
        <div className="lockup">
          <Logo size={56} />
          <h1>Create your account</h1>
          <p>Free — your sequences stay private to you.</p>
        </div>

        <div className="card">
          <div className="card-body">
            <form onSubmit={onSubmit}>
              <div className="field">
                <label htmlFor="name">Full name</label>
                <input
                  id="name"
                  type="text"
                  autoComplete="name"
                  value={name}
                  onChange={(e) => setName(e.target.value)}
                  placeholder="Dr Mohamed Merzoug"
                  required
                />
              </div>
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
                  autoComplete="new-password"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                  placeholder="At least 8 characters, letters and numbers"
                  required
                />
              </div>
              <div className="field">
                <label htmlFor="password2">Confirm password</label>
                <input
                  id="password2"
                  type="password"
                  autoComplete="new-password"
                  value={password2}
                  onChange={(e) => setPassword2(e.target.value)}
                  required
                />
              </div>

              {/* Carries the password rules when one is rejected, so it is
                  the text a reader most needs announced rather than found. */}
              {error && <div className="notice notice-error" role="alert">{error}</div>}

              <button className="btn btn-primary" type="submit" disabled={busy}>
                {busy && <span className="spinner" />}
                {busy ? "Creating…" : "Create account"}
              </button>
            </form>
          </div>
        </div>

        <p className="alt">
          Already registered? <Link to="/login">Sign in</Link>
        </p>
      </div>
    </div>
  );
}
