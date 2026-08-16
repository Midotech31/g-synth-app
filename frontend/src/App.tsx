import { Suspense, lazy } from "react";
import {
  BrowserRouter, Navigate, NavLink, Outlet, Route, Routes, useLocation,
} from "react-router-dom";

import { AuthProvider, useAuth } from "./auth/AuthContext";
import ErrorBoundary from "./components/ErrorBoundary";
import Icon, { type IconName } from "./components/Icon";
import { Logo } from "./components/Logo";
import Home from "./pages/Home";
import Login from "./pages/Login";
import Signup from "./pages/Signup";
import { WorkspaceStateProvider } from "./state/WorkspaceStateContext";

/**
 * The workspaces load on demand.
 *
 * Eagerly importing all of them put the sequence viewer and the trace
 * renderer — the two heaviest things here — into the bundle that every
 * signed-in user waits for, including on the pages that draw neither. Home,
 * Login and Signup stay eager: they are the first thing anyone sees, and
 * splitting them would trade a smaller bundle for a spinner on arrival.
 */
const Align = lazy(() => import("./pages/Align"));
const Clone = lazy(() => import("./pages/Clone"));
const Dashboard = lazy(() => import("./pages/Dashboard"));
const Design = lazy(() => import("./pages/Design"));
const Help = lazy(() => import("./pages/Help"));
const Learn = lazy(() => import("./pages/Learn"));
const Optimise = lazy(() => import("./pages/Optimise"));
const Pcr = lazy(() => import("./pages/Pcr"));
const Verify = lazy(() => import("./pages/Verify"));
const Viewer = lazy(() => import("./pages/Viewer"));

function Rail() {
  const { user, signOut } = useAuth();
  const initials = (user?.name || user?.email || "GS")
    .split(/\s+|@/)
    .filter(Boolean)
    .slice(0, 2)
    .map((part) => part[0]?.toUpperCase())
    .join("");
  const groups: Array<{
    label: string;
    items: Array<{ to: string; label: string; icon: IconName }>;
  }> = [
    {
      label: "Design",
      items: [
        { to: "/design", label: "Design", icon: "helix" },
        { to: "/optimise", label: "Optimise", icon: "target" },
      ],
    },
    {
      label: "Build",
      items: [
        { to: "/pcr", label: "PCR", icon: "microscope" },
        { to: "/clone", label: "Clone", icon: "plate" },
      ],
    },
    {
      label: "Verify",
      items: [
        { to: "/verify", label: "Check", icon: "check" },
        { to: "/align", label: "Compare", icon: "scales" },
        { to: "/learn", label: "Learn", icon: "book" },
      ],
    },
  ];

  return (
    <aside className="rail">
      <NavLink to="/" end className="brand" aria-label="G-Synth home">
        <Logo size={30} tagline={false} />
        <span className="ver">v3</span>
      </NavLink>

      <nav aria-label="Workspace">
        {groups.map((group) => (
          <div className="rail-group" key={group.label}>
            <div className="rail-group-label">{group.label}</div>
            {group.items.map((item) => (
              <NavLink
                key={item.to}
                to={item.to}
                className={({ isActive }) => (isActive ? "active" : "")}
              >
                <Icon name={item.icon} size={20} />
                <span>{item.label}</span>
              </NavLink>
            ))}
          </div>
        ))}
        <div className="rail-utilities">
          <NavLink to="/projects" className={({ isActive }) => (isActive ? "active" : "")}>
            <Icon name="book" size={20} />
            <span>Projects</span>
          </NavLink>
          <NavLink to="/help" className={({ isActive }) => (isActive ? "active" : "")}>
            <Icon name="target" size={20} />
            <span>Help</span>
          </NavLink>
        </div>
      </nav>

      <div className="spacer" />

      <div className="account">
        <div className="account-avatar" aria-hidden="true">{initials || "GS"}</div>
        <div className="who">
          <div className="n">{user?.name || "Signed in"}</div>
          <div className="e">{user?.email}</div>
        </div>
        <button className="rail-signout" onClick={() => void signOut()} title="Sign out">
          <span className="sr-only">Sign out</span>
          <Icon name="arrowRight" size={18} />
        </button>
      </div>
    </aside>
  );
}

/** Renders the app shell, or bounces to /login when signed out. */
function Protected() {
  const { user, loading } = useAuth();
  // Read before the early returns below: a hook called after a conditional
  // return runs on some renders and not others, which is exactly the
  // ordering React forbids.
  const { pathname } = useLocation();

  if (loading) {
    return (
      <div className="center-note" role="status" aria-live="polite">
        <span className="spinner" />
        <span>Loading…</span>
      </div>
    );
  }
  if (!user) return <Navigate to="/login" replace />;

  return (
    <WorkspaceStateProvider>
      <div className="shell">
      {/* First thing in the tab order, so the nine rail links can be passed
          over. `tabIndex` on the target because following a fragment moves
          the caret but not focus in several browsers. */}
      <a className="skip-link" href="#main">
        Skip to main content
      </a>
      <Rail />
      <main className="canvas" id="main" role="main" tabIndex={-1}>
        {/* Inside the shell, so a page that fails to render leaves the rail
            standing and the reader can move somewhere else. Keyed on the
            path: without that, the boundary stays latched after a
            navigation and the next page renders as the previous error. */}
        <ErrorBoundary key={pathname}>
          <Suspense
            fallback={
              <div className="center-note" role="status" aria-live="polite">
                <span className="spinner" />
                <span>Loading…</span>
              </div>
            }
          >
            <Outlet />
          </Suspense>
        </ErrorBoundary>
      </main>
      </div>
    </WorkspaceStateProvider>
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
            <Route path="/pcr" element={<Pcr />} />
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

