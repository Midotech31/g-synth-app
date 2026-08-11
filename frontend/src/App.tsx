import { BrowserRouter, Navigate, NavLink, Outlet, Route, Routes } from "react-router-dom";

import { AuthProvider, useAuth } from "./auth/AuthContext";
import Align from "./pages/Align";
import Clone from "./pages/Clone";
import Dashboard from "./pages/Dashboard";
import Design from "./pages/Design";
import Help from "./pages/Help";
import Home from "./pages/Home";
import Learn from "./pages/Learn";
import Optimise from "./pages/Optimise";
import Verify from "./pages/Verify";
import Login from "./pages/Login";
import Signup from "./pages/Signup";
import Viewer from "./pages/Viewer";
import { Logo } from "./components/Logo";

function Rail() {
  const { user, signOut } = useAuth();
  return (
    <aside className="rail">
      <div className="brand">
        <Logo size={30} tagline={false} />
        <span className="ver">v3</span>
      </div>

      <nav>
        <NavLink to="/" end className={({ isActive }) => (isActive ? "active" : "")}>
          Home
        </NavLink>
        <NavLink to="/design" className={({ isActive }) => (isActive ? "active" : "")}>
          Design
        </NavLink>
        <NavLink to="/optimise" className={({ isActive }) => (isActive ? "active" : "")}>
          Optimise
        </NavLink>
        <NavLink to="/clone" className={({ isActive }) => (isActive ? "active" : "")}>
          Clone
        </NavLink>
        <NavLink to="/verify" className={({ isActive }) => (isActive ? "active" : "")}>
          Check
        </NavLink>
        <NavLink to="/align" className={({ isActive }) => (isActive ? "active" : "")}>
          Compare
        </NavLink>
        <NavLink to="/learn" className={({ isActive }) => (isActive ? "active" : "")}>
          Learn
        </NavLink>
        <NavLink to="/projects" className={({ isActive }) => (isActive ? "active" : "")}>
          Projects
        </NavLink>
        <NavLink to="/help" className={({ isActive }) => (isActive ? "active" : "")}>
          Help
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
            <Route path="/" element={<Home />} />
            <Route path="/design" element={<Design />} />
            <Route path="/optimise" element={<Optimise />} />
            <Route path="/clone" element={<Clone />} />
            <Route path="/verify" element={<Verify />} />
            <Route path="/align" element={<Align />} />
            <Route path="/learn" element={<Learn />} />
            <Route path="/projects" element={<Dashboard />} />
            <Route path="/projects/:id" element={<Viewer />} />
            <Route path="/help" element={<Help />} />
          </Route>
          <Route path="*" element={<Navigate to="/" replace />} />
        </Routes>
      </AuthProvider>
    </BrowserRouter>
  );
}
