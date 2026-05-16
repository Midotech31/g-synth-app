# G-Synth 3.0 — production Docker image
FROM python:3.12-slim AS base

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# System dependencies — minimal: gcc only for primer3-py wheel fallback
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Optional binaries for MSA (uncomment to enable):
# RUN apt-get update && apt-get install -y --no-install-recommends \
#        mafft muscle clustalo \
#    && rm -rf /var/lib/apt/lists/*

# Copy project metadata AND source before install. setuptools needs the
# `gsynth_core` / `gsynth_ui` packages on disk to discover them via
# pyproject.toml's `packages.find` configuration; installing before COPY
# would silently install only the dependencies and leave the package empty.
COPY pyproject.toml README.md app.py ./
COPY assets/      ./assets/
COPY gsynth_core/ ./gsynth_core/
COPY gsynth_ui/   ./gsynth_ui/
COPY pages/       ./pages/
COPY .streamlit/  ./.streamlit/

RUN pip install --upgrade pip && \
    pip install ".[ui]"

# Non-root user
RUN useradd --create-home --shell /bin/bash gsynth && \
    chown -R gsynth:gsynth /app
USER gsynth

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --start-period=20s --retries=3 \
    CMD python -c "import urllib.request, sys; \
        sys.exit(0 if urllib.request.urlopen('http://localhost:8501/_stcore/health').status == 200 else 1)"

ENTRYPOINT ["streamlit", "run", "app.py", \
            "--server.port=8501", "--server.address=0.0.0.0"]
