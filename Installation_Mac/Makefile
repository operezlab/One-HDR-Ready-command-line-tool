# ============================================================
# One-HDR-Ready + CRISPOR Installer (Python 3.9 venv)
# Cross-platform: Ubuntu/WSL + macOS
# ============================================================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c

PYTHON := python3.9
VENV_DIR := venv
ACTIVATE := . $(VENV_DIR)/bin/activate

OS := $(shell uname)

# Default target
all: install

# ------------------------------------------------------------
# Step 1. Check Python 3.9 installed (Ubuntu + macOS)
# ------------------------------------------------------------
check_python:
	@echo "Checking for Python 3.9..."
	if ! command -v $(PYTHON) >/dev/null 2>&1; then \
		echo "Python 3.9 not found."; \
		if [ "$(OS)" = "Linux" ]; then \
			echo "Installing Python 3.9 for Ubuntu/WSL (requires sudo)..."; \
			sudo apt update && sudo apt install -y python3.9 python3.9-venv; \
		elif [ "$(OS)" = "Darwin" ]; then \
			echo "Installing Python 3.9 for macOS via Homebrew..."; \
			brew install python@3.9; \
		fi; \
	else \
		echo "Python 3.9 found."; \
	fi

# ------------------------------------------------------------
# Step 2. Create Python 3.9 venv
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
# Step 3. Install CRISPOR-safe dependencies (NO TF, NO Keras)
# ------------------------------------------------------------
deps: venv
	@echo "Installing pinned CRISPOR dependencies..."
	$(ACTIVATE); \
	pip install \
		biopython==1.79 \
		numpy==1.23.5 \
		pandas==1.5.3 \
		scikit-learn==1.0.2 \
		matplotlib==3.7.1 \
		twobitreader==3.1.7 \
		xlwt==1.3.0 \
		h5py==3.8.0 \
		pytabix==0.1 \
		rs3==0.0.18 \
		lmdbm==0.0.6

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
	@echo "CRISPOR installed at ./crisporWebsite"

# ------------------------------------------------------------
# Step 5. Download hg38 genome index
# ------------------------------------------------------------
hg38:
	@echo "Downloading hg38 genome index..."
	mkdir -p crisporGenomes/hg38
	cd crisporGenomes/hg38; \
	for f in hg38.2bit hg38.fa hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed genomeInfo.tab; do \
		if [ ! -f $$f ]; then \
			echo "Downloading $$f..."; \
			curl -sO https://crispor.gi.ucsc.edu/genomes/hg38/$$f; \
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
	@echo "Scripts ready in ./One-HDR-Ready-command-line-tool"

# ------------------------------------------------------------
# Step 7. Full installation pipeline
# ------------------------------------------------------------
install: check_python venv deps crispor hg38 onehdr
	@echo "==========================================================="
	@echo " Installation complete!"
	@echo ""
	@echo " To activate the environment:"
	@echo "     source $(VENV_DIR)/bin/activate"
	@echo ""
	@echo " To run One-HDR-Ready:"
	@echo "     python One-HDR-Ready-command-line-tool/OneHDRReady.py"
	@echo "==========================================================="

# ------------------------------------------------------------
# Optional cleanup
# ------------------------------------------------------------
clean:
	@echo "Removing local installs..."
	rm -rf $(VENV_DIR) crisporWebsite crisporGenomes One-HDR-Ready-command-line-tool
