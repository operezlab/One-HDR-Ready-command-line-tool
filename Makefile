# ============================================================
# One-HDR-Ready Installer
# - Python deps
# - CRISPOR
# - hg38 genome
# - Azimuth retraining (safe + idempotent)
#
# NOTE:
#   Conda environment MUST already be activated
# ============================================================

REQ_FILE     = requirements.txt
CRISPOR_DIR  = crisporWebsite
GENOME_DIR   = crisporGenomes
AZI_DIR      = $(CRISPOR_DIR)/bin/Azimuth-2.0
AZI_MARKER   = $(AZI_DIR)/.trained
PYTHON       = python

# Default target
all: install

# ------------------------------------------------------------
# Install Python dependencies
# ------------------------------------------------------------
install_requirements:
	@echo ">>> Installing Python dependencies into active conda environment"
	pip install -r $(REQ_FILE)

	@echo ">>> Installing NCBI BLAST+"
	conda install -y -c bioconda blast

# ------------------------------------------------------------
# Install CRISPOR
# ------------------------------------------------------------
crispor:
	@echo ">>> Installing CRISPOR"
	if [ ! -d "$(CRISPOR_DIR)" ]; then \
		git clone https://github.com/maximilianh/crisporWebsite.git $(CRISPOR_DIR); \
	else \
		echo "CRISPOR already installed — skipping clone."; \
	fi

# ------------------------------------------------------------
# Install hg38 genome
# ------------------------------------------------------------
hg38:
	@echo ">>> Installing hg38 genome index"
	mkdir -p $(GENOME_DIR)/hg38
	cd $(GENOME_DIR)/hg38 && \
	for f in hg38.2bit hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed; do \
		if [ ! -f $$f ]; then \
			echo "Downloading $$f"; \
			curl -s -O https://crispor.gi.ucsc.edu/genomes/hg38/$$f; \
		else \
			echo "$$f already exists"; \
		fi; \
	done

# ------------------------------------------------------------
# Retrain Azimuth models (required for sklearn/numpy compatibility)
# ------------------------------------------------------------
azimuth-retrain:
	@if [ -f "$(AZI_MARKER)" ]; then \
		echo ">>> Azimuth already trained — skipping."; \
	else \
		echo "======================================================="; \
		echo " Retraining Azimuth models for current environment"; \
		echo "======================================================="; \
		if [ ! -d "$(AZI_DIR)" ]; then \
			echo "ERROR: Azimuth directory not found at $(AZI_DIR)"; \
			exit 1; \
		fi; \
		cd $(AZI_DIR) && $(PYTHON) model_comparison.py && touch $(AZI_MARKER); \
		echo ">>> Azimuth retraining complete"; \
	fi

# ------------------------------------------------------------
# Full installation
# ------------------------------------------------------------
install: install_requirements crispor hg38 azimuth-retrain
	@echo "======================================================="
	@echo " One-HDR-Ready installation complete"
	@echo " Active conda environment: $$CONDA_DEFAULT_ENV"
	@echo " CRISPOR directory: $(CRISPOR_DIR)"
	@echo " Genome directory:  $(GENOME_DIR)"
	@echo "======================================================="

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------
clean:
	rm -rf $(CRISPOR_DIR) $(GENOME_DIR)
