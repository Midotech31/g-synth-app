"""Test configuration. Shims pytest if it isn't installed."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

try:
    import pytest  # noqa: F401
except ImportError:
    import types

    class _Raises:
        def __init__(self, exc): self.exc = exc
        def __enter__(self): return self
        def __exit__(self, et, ev, tb):
            if et is None:
                raise AssertionError(f"Did not raise {self.exc}")
            return issubclass(et, self.exc)

    pytest = types.ModuleType("pytest")
    pytest.raises = lambda exc: _Raises(exc)
    pytest.fixture = lambda *a, **kw: (a[0] if a and callable(a[0]) else (lambda f: f))
    pytest.skip = lambda *a, **kw: None
    class _Mark:
        skip = staticmethod(lambda *a, **kw: (lambda fn: fn))
        parametrize = staticmethod(lambda *a, **kw: (lambda fn: fn))
    pytest.mark = _Mark()
    sys.modules["pytest"] = pytest
