import { BrowserRouter, Navigate, NavLink, Outlet, Route, Routes } from "react-router-dom";

import { AuthProvider, useAuth } from "./auth/AuthContext";
import Dashboard from "./pages/Dashboard";
import Login from "./pages/Login";
import Signup from "./pages/Signup";
import Viewer from "./pages/Viewer";

function Rail() {
  const { user, signOut } = useAuth();
  return (
    <aside className="rail">
      <div className="brand">
        <span className="mark">🧬</span>
        <span>
          <span className="name">G-Synth</span>
          <br />
          <span className="ver">v3 · workspace</span>
        </span>
      </div>

      <nav>
        <NavLink to="/" end className={({ isActive }) => (isActive ? "active" : "")}>
          Projects
        </NavLink>
      </nav>

      <div className="spacer" />

      <div className="account">
        <div className="who">
          <div className="n">{user?.name || "Signed in"}</div>
          <div className="e">{user?.email}</div>
        </div>
        <button className="btn btn-rail" onClick={() => void signOut()}>
          Sign out
        </button>
      </div>
    </aside>
  );
}

/** Renders the app shell, or bounces to /login when signed out. */
function Protected() {
  const { user, loading } = useAuth();

  if (loading) {
    return (
      <div className="center-note">
        <span className="spinner" />
        <span>Loading…</span>
      </div>
    );
  }
  if (!user) return <Navigate to="/login" replace />;

  return (
    <div className="shell">
      <Rail />
      <main className="canvas">
        <Outlet />
      </main>
    </div>
  );
}

/** Keeps signed-in users away from the auth screens. */
function PublicOnly({ children }: { children: React.ReactNode }) {
  const { user, loading } = useAuth();
  if (loading) return null;
  if (user) return <Navigate to="/" replace />;
  return <>{children}</>;
}

export default function App() {
  return (
    <BrowserRouter>
      <AuthProvider>
        <Routes>
          <Route
            path="/login"
            element={
              <PublicOnly>
                <Login />
              </PublicOnly>
            }
          />
          <Route
            path="/signup"
            element={
              <PublicOnly>
                <Signup />
              </PublicOnly>
            }
          />
          <Route element={<Protected />}>
            <Route path="/" element={<Dashboard />} />
            <Route path="/projects/:id" element={<Viewer />} />
          </Route>
          <Route path="*" element={<Navigate to="/" replace />} />
        </Routes>
      </AuthProvider>
    </BrowserRouter>
  );
}
