# =============================================================================
# Active Geometry: Reproducible Verification Container
# =============================================================================
# Build:  docker build -t active-geometry .
# Run:    docker run --rm active-geometry
# =============================================================================

FROM python:3.11-slim

# Environment
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    ELAN_HOME="/root/.elan" \
    PATH="/root/.elan/bin:$PATH"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash ca-certificates curl git jq libgmp-dev \
  && rm -rf /var/lib/apt/lists/*

# Install Lean 4 via elan
RUN curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | \
    sh -s -- -y --default-toolchain leanprover/lean4:stable

WORKDIR /repo

# Python dependencies (cached layer)
COPY requirements.txt /repo/requirements.txt
RUN pip install --no-cache-dir -r /repo/requirements.txt && \
    pip install --no-cache-dir jupyter nbconvert nbclient nbformat sympy mpmath

# Copy repository
COPY . /repo

# Make scripts executable
RUN chmod +x /repo/run_all_verifications.sh 2>/dev/null || true

# Artifacts mountpoint
VOLUME ["/artifacts"]

# Entrypoint
ENTRYPOINT ["/bin/bash", "-lc", "cd /repo && ./run_all_verifications.sh"]
