# ============================================================
# One-HDR-Ready Web Tool Installer (Python 3.9 venv, local)
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
# Step 1. Ensure Python 3.9 and venv are installed
# ------------------------------------------------------------
check_python:
	@echo "Checking for Python 3.9..."
	if ! command -v $(PYTHON) >/dev/null 2>&1; then \
		echo "Python 3.9 not found. Installing (requires sudo)..."; \
		sudo apt update && sudo apt install -y python3.9 python3.9-venv; \
	else \
		echo "Python 3.9 found."; \
	fi

# ------------------------------------------------------------
# Step 2. Create and activate virtual environment
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
# Step 3. Install Python dependencies
# ------------------------------------------------------------
deps: venv
	@echo "Installing Python dependencies..."
	$(ACTIVATE); pip install biopython pandas numpy flask gunicorn requests matplotlib primer3-py gitpython wget

# ------------------------------------------------------------
# Step 4. Install CRISPOR (CLI only)
# ------------------------------------------------------------
crispor: deps
	@echo "Installing CRISPOR..."
	if [ ! -d "crisporWebsite" ]; then \
		git clone https://github.com/maximilianh/crisporWebsite.git crisporWebsite; \
	else \
		echo "CRISPOR already cloned."; \
	fi
	$(ACTIVATE); pip install biopython numpy scikit-learn pandas twobitreader xlwt keras tensorflow h5py rs3 pytabix lmdbm
	@echo "CRISPOR CLI ready at ./crisporWebsite/crispor.py."

# ------------------------------------------------------------
# Step 5. Download hg38 genome index
# ------------------------------------------------------------
hg38:
	@echo "Downloading hg38 genome index..."
	mkdir -p crisporGenomes/hg38
	cd crisporGenomes/hg38
	for f in hg38.2bit hg38.fa hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed genomeInfo.tab; do \
		if [ ! -f $$f ]; then \
			echo "Downloading $$f..."; \
			wget -q https://crispor.gi.ucsc.edu/genomes/hg38/$$f; \
		else \
			echo "$$f already exists."; \
		fi; \
	done

# ------------------------------------------------------------
# Step 6. Clone One-HDR-Ready tool
# ------------------------------------------------------------
onehdr:
	@echo "Cloning One-HDR-Ready command-line tool..."
	if [ ! -d "One-HDR-Ready-command-line-tool" ]; then \
		git clone https://github.com/operezlab/One-HDR-Ready-command-line-tool.git One-HDR-Ready-command-line-tool; \
	else \
		echo "Repository already cloned."; \
	fi
	@echo "Scripts ready in ./One-HDR-Ready-command-line-tool."

# ------------------------------------------------------------
# Step 7. Aggregate installation
# ------------------------------------------------------------
install: check_python venv deps crispor hg38 onehdr
	@echo "==========================================================="
	@echo "Installation complete."
	@echo "To use the environment, run:"
	@echo "  source $(VENV_DIR)/bin/activate"
	@echo "Then launch the tool with:"
	@echo "  python One-HDR-Ready-command-line-tool/OneHDRReady.py"
	@echo "==========================================================="

# ------------------------------------------------------------
# Optional cleanup
# ------------------------------------------------------------
clean:
	@echo "Removing local installations..."
	rm -rf $(VENV_DIR) crisporWebsite crisporGenomes One-HDR-Ready-command-line-tool
