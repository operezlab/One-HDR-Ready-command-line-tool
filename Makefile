# ============================================================
# One-HDR-Ready Web Tool Full Installer (Ubuntu / WSL)
# ============================================================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -e -o pipefail -c

# Default target
all: install

# ------------------------------------------------------------
# Step 1. Install Miniconda if missing
# ------------------------------------------------------------
install_miniconda:
	@echo "🔍 Checking Miniconda..."
	if ! command -v conda >/dev/null 2>&1; then
		echo "📦 Installing Miniconda..."
		wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/Miniconda3-latest-Linux-x86_64.sh
		bash ~/Miniconda3-latest-Linux-x86_64.sh -b -p $$HOME/miniconda3
		echo "✅ Miniconda installed."
	else
		echo "✅ Miniconda already installed."
	fi

# ------------------------------------------------------------
# Step 2. Create and activate environment + dependencies
# ------------------------------------------------------------
env:
	@echo "🧱 Setting up conda environment 'crispr-webtool'..."
	export PATH="$$HOME/miniconda3/bin:$$PATH"
	eval "$$(conda shell.bash hook)"
	if conda env list | grep -q "crispr-webtool"; then
		echo "🔁 Environment already exists. Skipping creation."
	else
		conda create -y -n crispr-webtool python=3.9
	fi
	conda activate crispr-webtool
	conda install -y biopython pandas numpy flask gunicorn requests matplotlib git wget curl make
	pip install primer3-py

# ------------------------------------------------------------
# Step 3. Install CRISPOR
# ------------------------------------------------------------
crispor:
	@echo "Installing CRISPOR..."
	export PATH="$$HOME/miniconda3/bin:$$PATH"
	eval "$$(conda shell.bash hook)"
	conda activate crispr-webtool
	if [ ! -d "$$HOME/crisporWebsite" ]; then
		git clone https://github.com/maximilianh/crisporWebsite.git $$HOME/crisporWebsite
	else
		echo "CRISPOR already cloned."
	fi
	@echo "Installing CRISPOR Python dependencies..."
	pip install biopython numpy scikit-learn pandas twobitreader xlwt keras tensorflow h5py rs3 pytabix lmdbm
	@echo "CRISPOR CLI setup complete. Skipping 'make' step (web interface not required)."

# ------------------------------------------------------------
# Step 4. Download hg38 genome index
# ------------------------------------------------------------
hg38:
	@echo "📥 Downloading hg38 genome index..."
	mkdir -p $$HOME/crisporGenomes/hg38
	cd $$HOME/crisporGenomes/hg38
	for f in hg38.2bit hg38.fa hg38.fa.amb hg38.fa.ann hg38.fa.bwt hg38.fa.pac hg38.fa.sa hg38.sizes hg38.segments.bed genomeInfo.tab; do
		if [ ! -f $$f ]; then
			echo "Fetching $$f..."
			wget -q https://crispor.gi.ucsc.edu/genomes/hg38/$$f
		else
			echo "$$f already exists."
		fi
	done

# ------------------------------------------------------------
# Step 5. Clone One-HDR-Ready tool
# ------------------------------------------------------------
onehdr:
	@echo "Cloning One-HDR-Ready command-line tool..."
	export PATH="$$HOME/miniconda3/bin:$$PATH"
	eval "$$(conda shell.bash hook)"
	conda activate crispr-webtool
	if [ ! -d "$$HOME/One-HDR-Ready-command-line-tool" ]; then
		git clone https://github.com/operezlab/One-HDR-Ready-command-line-tool.git $$HOME/One-HDR-Ready-command-line-tool
	else
		echo "Repository already cloned."
	fi
	@echo "One-HDR-Ready cloned successfully. Skipping pip installation (no setup.py found)."

# ------------------------------------------------------------
# Step 6. Aggregate target
# ------------------------------------------------------------
install: install_miniconda env crispor hg38 onehdr
	@echo "==========================================================="
	@echo "🎉 Installation complete!"
	@echo "To launch the One-HDR-Ready Web Tool, run:"
	@echo "  conda activate crispr-webtool"
	@echo "  python3 OneHDRReady.py"
	@echo "==========================================================="

# ------------------------------------------------------------
# Optional cleanup
# ------------------------------------------------------------
clean:
	@echo "🧹 Cleaning installations..."
	rm -rf $$HOME/One-HDR-Ready-command-line-tool $$HOME/crisporWebsite $$HOME/crisporGenomes
