# ============================================================
# One-HDR-Ready Web Tool Installer (macOS, Python 3.9 venv)
# ============================================================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c

PYTHON := python3.9
VENV_DIR := venv
ACTIVATE := . $(VENV_DIR)/bin/activate

# Default target
all: install

# ------------------------------------------------------------
# Step 1 — Ensure Python 3.9 is installed via Homebrew
# ------------------------------------------------------------
check_python:
	@echo "Checking for Python 3.9..."
	if ! command -v $(PYTHON) >/dev/null 2>&1; then \
		echo "Python 3.9 not found. Installing with Homebrew..."; \
		if ! command -v brew >/dev/null 2>&1; then \
			echo "Homebrew not found. Install it from https://brew.sh"; \
			exit 1; \
		fi; \
		brew update && brew install python@3.9 && brew link python@3.9 --force --overwrite; \
	else \
		echo "Python 3.9 found."; \
	fi

# ------------------------------------------------------------
# Step 2 — Create and activate virtual environment
# ------------------------------------------------------------
venv: check_python
	@echo "Creating Python 3.9 virtual environment..."
	if [ ! -d "$(VENV_DIR)" ]; then \
		$(PYTHON) -m venv $(VENV_DIR); \
		echo "Virtual environment created at $(VENV_DIR)."; \
	else \
		echo "Virtual environment already exists."; \
	fi
	$(ACTIVATE); pip install --upgrade pip setuptools wheel

# ------------------------------------------------------------
# Step 3 — Install Python dependencies (NO LightGBM)
# ------------------------------------------------------------
deps: venv
	@echo "Installing required Python packages..."
	$(ACTIVATE); \
	pip install biopython numpy scikit-learn pandas twobitreader xlwt keras tensorflow h5py pytabix matplotlib lmdbm gitpython; \
	pip install rs3 --no-deps

# ------------------------------------------------------------
# Step 4 — Install CRISPOR
# ------------------------------------------------------------
crispor: deps
	@echo "Installing CRISPOR..."
	if [ ! -d "crisporWebsite" ]; then \
		git clone https://github.com/maximilianh/crisporWebsite.git crisporWebsite; \
	else \
		echo "CRISPOR already cloned."; \
	fi
	@echo "CRISPOR CLI ready at ./crisporWebsite/crispor.py."

# ------------------------------------------------------------
# Step 5 — Download hg38 genome index (curl version)
# ------------------------------------------------------------
hg38:
	@echo "Downloading hg38 genome index..."
	mkdir -p crisporGenomes/hg38
	cd crisporGenomes/hg38; \
	for f in hg38.2bit hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed genomeInfo.tab; do \
		if [ ! -f $$f ]; then \
			echo "Downloading $$f..."; \
			curl -fsSLO https://crispor.gi.ucsc.edu/genomes/hg38/$$f; \
		else \
			echo "$$f already exists."; \
		fi; \
	done

# ------------------------------------------------------------
# Step 6 — Clone One-HDR-Ready tool
# ------------------------------------------------------------
onehdr:
	@echo "Cloning One-HDR-Ready command-line tool..."
	if [ ! -d "One-HDR-Ready-command-line-tool" ]; then \
		git clone https://github.com/operezlab/One-HDR-Ready-command-line-tool.git One-HDR-Ready-command-line-tool; \
	else \
		echo "Repository already cloned."; \
	fi
	@echo "Scripts available in ./One-HDR-Ready-command-line-tool."

# ------------------------------------------------------------
# Step 7 — Aggregate installation
# ------------------------------------------------------------
install: check_python venv deps crispor hg38 onehdr
	@echo "==========================================================="
	@echo "Installation complete!"
	@echo ""
	@echo "To activate the environment, run:"
	@echo "  source $(VENV_DIR)/bin/activate"
	@echo ""
	@echo "To launch the One-HDR-Ready tool, run:"
	@echo "  python One-HDR-Ready-command-line-tool/OneHDRReady.py"
	@echo "==========================================================="

# ------------------------------------------------------------
# Optional cleanup
# ------------------------------------------------------------
clean:
	@echo "Removing local installations..."
	rm -rf $(VENV_DIR) crisporWebsite crisporGenomes One-HDR-Ready-command-line-tool
