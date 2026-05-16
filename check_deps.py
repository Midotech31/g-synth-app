# check_deps.py
import importlib
import sys

# read your requirements.txt
with open("requirements.txt") as f:
    reqs = [line.strip().split("==")[0].split(">=")[0] for line in f if line.strip() and not line.startswith("#")]

missing = []
for pkg in reqs:
    try:
        importlib.import_module(pkg)
    except ImportError:
        missing.append(pkg)

if missing:
    print("❌ Missing or failed imports for:", ", ".join(missing))
    sys.exit(1)
else:
    print("✅ All required packages imported successfully!")
    sys.exit(0)