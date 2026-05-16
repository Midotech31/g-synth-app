#!/usr/bin/env python3
"""
Simple test runner for environments without pytest installed.

Discovers `test_*.py` files in `tests/` and runs every top-level callable
named `test_*` plus methods named `test_*` on classes named `Test*`.

Prefer `pytest` when available (`pip install pytest`).
"""
from __future__ import annotations
import importlib.util
import inspect
import sys
import traceback
from pathlib import Path


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

# Load conftest first so the pytest shim is in place
_conftest = ROOT / "tests" / "conftest.py"
if _conftest.exists():
    import importlib.util as _il
    _spec = _il.spec_from_file_location("conftest", _conftest)
    _mod = _il.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)


def _load_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _discover(test_dir: Path) -> list[tuple[str, callable]]:
    tests: list[tuple[str, callable]] = []
    for test_file in sorted(test_dir.glob("test_*.py")):
        try:
            mod = _load_module(test_file)
        except Exception as e:
            tests.append((f"{test_file.stem}::IMPORT_ERROR",
                          lambda e=e: (_ for _ in ()).throw(e)))
            continue
        for name, obj in vars(mod).items():
            if name.startswith("test_") and callable(obj):
                tests.append((f"{test_file.stem}::{name}", obj))
            elif name.startswith("Test") and inspect.isclass(obj):
                for mname, meth in inspect.getmembers(obj, inspect.isfunction):
                    if mname.startswith("test_"):
                        instance = obj()
                        tests.append((f"{test_file.stem}::{name}::{mname}",
                                      getattr(instance, mname)))
    return tests


def main() -> int:
    test_dir = ROOT / "tests"
    if not test_dir.exists():
        print("No tests/ directory found")
        return 1
    tests = _discover(test_dir)
    if not tests:
        print("No tests collected")
        return 1
    passed = failed = errored = 0
    failures: list[tuple[str, str]] = []
    for name, fn in tests:
        try:
            fn()
            passed += 1
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            failures.append((name, traceback.format_exc()))
            print(f"  FAIL  {name}: {e}")
        except Exception as e:
            errored += 1
            failures.append((name, traceback.format_exc()))
            print(f"  ERROR {name}: {type(e).__name__}: {e}")
    total = len(tests)
    print()
    print(f"==== {total} tests:  {passed} passed, {failed} failed, {errored} errors ====")
    if failures:
        print()
        print("DETAILS:")
        for name, tb in failures:
            print(f"\n--- {name} ---\n{tb}")
    return 0 if (failed == 0 and errored == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
