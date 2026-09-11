import { useState } from "react";
import { Link, useNavigate } from "react-router-dom";

import { ApiError } from "../api/client";
import { useAuth } from "../auth/AuthContext";
import AuthShell from "../components/AuthShell";
import PasswordInput from "../components/PasswordInput";

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
    <AuthShell
      eyebrow="Private scientific workspace"
      title="Create your account"
      description="Start a private workspace for synthesis-ready designs and their validation evidence."
      footer={<>Already registered? <Link to="/login">Sign in</Link></>}
    >
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
                <PasswordInput
                  id="password"
                  label="Password"
                  autoComplete="new-password"
                  value={password}
                  onChange={setPassword}
                  placeholder="At least 8 characters, letters and numbers"
                />
              </div>
              <div className="field">
                <label htmlFor="password2">Confirm password</label>
                <PasswordInput
                  id="password2"
                  label="Confirm password"
                  autoComplete="new-password"
                  value={password2}
                  onChange={setPassword2}
                />
              </div>


              {error && <div className="notice notice-error" role="alert">{error}</div>}

              <button className="btn btn-primary" type="submit" disabled={busy}>
                {busy && <span className="spinner" />}
                {busy ? "Creating…" : "Create account"}
              </button>
      </form>
    </AuthShell>
  );
}
