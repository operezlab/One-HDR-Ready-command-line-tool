# ============================================================
# One-HDR-Ready Installer (Dependencies + CRISPOR + hg38)
# Environment must already be activated!
# ============================================================

REQ_FILE = requirements.txt
CRISPOR_DIR = crisporWebsite
GENOME_DIR = crisporGenomes

# Default target
all: install

# ------------------------------------------------------------
# Install Python dependencies into current environment
# ------------------------------------------------------------
install_requirements:
	@echo ">>> Installing Python dependencies into active conda environment"
	pip install -r $(REQ_FILE)

	@echo "Installing NCBI BLAST+ ..."
	conda install -y -c bioconda blast

# ------------------------------------------------------------
# Install CRISPOR Website
# ------------------------------------------------------------
crispor:
	@echo "Installing CRISPOR..."
	if [ ! -d "$(CRISPOR_DIR)" ]; then \
		git clone https://github.com/maximilianh/crisporWebsite.git $(CRISPOR_DIR); \
	else \
		echo "CRISPOR already cloned."; \
	fi

# ------------------------------------------------------------
# Install HG38 genome database using curl
# ------------------------------------------------------------
hg38:
	@echo "Downloading hg38 genome index..."
	mkdir -p $(GENOME_DIR)/hg38
	cd $(GENOME_DIR)/hg38 && \
	for f in hg38.2bit hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed; do \
		if [ ! -f $$f ]; then \
			echo "Downloading $$f..."; \
			curl -s -O https://crispor.gi.ucsc.edu/genomes/hg38/$$f; \
		else \
			echo "$$f already exists."; \
		fi; \
	done

# ------------------------------------------------------------
# Full installation
# ------------------------------------------------------------
install: install_requirements crispor hg38
	@echo "======================================================="
	@echo " Installation complete!"
	@echo " Active environment:"
	@echo "     $$CONDA_DEFAULT_ENV"
	@echo " CRISPOR installed in: $(CRISPOR_DIR)"
	@echo " Genome files in: $(GENOME_DIR)"
	@echo "======================================================="

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------
clean:
	rm -rf $(CRISPOR_DIR) $(GENOME_DIR)
