import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useRef,
  useState,
  type Dispatch,
  type ReactNode,
  type SetStateAction,
} from "react";

type WorkspaceStore = Record<string, unknown>;

type WorkspaceStateContextValue = {
  store: WorkspaceStore;
  setStore: Dispatch<SetStateAction<WorkspaceStore>>;
};

const WorkspaceStateContext = createContext<WorkspaceStateContextValue | null>(null);

/**
 * Keeps unfinished workspace work alive while authenticated routes swap in
 * and out. The provider deliberately lives inside the protected shell: a
 * sign-out destroys it, so sequences and results cannot cross accounts.
 */
export function WorkspaceStateProvider({ children, identity = "anonymous" }: { children: ReactNode; identity?: string }) {
  const storageKey = `gsynth.workspace.${identity}`;
  const [store, setStore] = useState<WorkspaceStore>(() => {
    try {
      const saved = window.sessionStorage.getItem(storageKey);
      return saved ? JSON.parse(saved) as WorkspaceStore : {};
    } catch {
      return {};
    }
  });

  useEffect(() => {
    // Recompute derived results after restoring editable inputs.
    const draftEntries = Object.entries(store).filter(([key]) => (
      !/\.(result|report|primers|ligation|saved|selected|traceFiles|project)$/.test(key)
    ));
    try {
      window.sessionStorage.setItem(storageKey, JSON.stringify(Object.fromEntries(draftEntries)));
    } catch {
      // Session persistence is optional.
    }
  }, [storageKey, store]);

  return (
    <WorkspaceStateContext.Provider value={{ store, setStore }}>
      {children}
    </WorkspaceStateContext.Provider>
  );
}

/** A useState-compatible value that survives navigation between workspaces. */
export function useWorkspaceState<T>(
  key: string,
  initialValue: T,
): [T, Dispatch<SetStateAction<T>>, () => void] {
  const context = useContext(WorkspaceStateContext);
  if (!context) {
    throw new Error("useWorkspaceState must be used inside WorkspaceStateProvider");
  }

  const { store, setStore } = context;
  const initialRef = useRef(initialValue);
  const hasValue = Object.prototype.hasOwnProperty.call(store, key);
  const value = (hasValue ? store[key] : initialRef.current) as T;

  const setValue = useCallback<Dispatch<SetStateAction<T>>>((next) => {
    setStore((currentStore) => {
      const current = (Object.prototype.hasOwnProperty.call(currentStore, key)
        ? currentStore[key]
        : initialRef.current) as T;
      const resolved = typeof next === "function"
        ? (next as (previous: T) => T)(current)
        : next;
      return { ...currentStore, [key]: resolved };
    });
  }, [key, setStore]);

  const clearValue = useCallback(() => {
    setStore((currentStore) => {
      if (!Object.prototype.hasOwnProperty.call(currentStore, key)) return currentStore;
      const next = { ...currentStore };
      delete next[key];
      return next;
    });
  }, [key, setStore]);

  return [value, setValue, clearValue];
}
