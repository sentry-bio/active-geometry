# =============================================================================
# Active Geometry: Master Makefile
# =============================================================================
# Targets for training, validation, and figure generation
# =============================================================================

.PHONY: help install test verify-lean train-all-seeds analyze-convergence \
        viral-validation tree-validation figures fetch-data clean validate-all

# Default target
help:
	@echo "Active Geometry - Makefile Targets"
	@echo "==================================="
	@echo ""
	@echo "Setup:"
	@echo "  make install           Install Python dependencies"
	@echo "  make test              Run smoke tests"
	@echo ""
	@echo "Theory Verification:"
	@echo "  make verify-lean       Compile Lean 4 proofs"
	@echo ""
	@echo "Training (requires GPU):"
	@echo "  make train-seed-0      Train with seed=0"
	@echo "  make train-seed-42     Train with seed=42"
	@echo "  make train-seed-123    Train with seed=123"
	@echo "  make train-seed-456    Train with seed=456"
	@echo "  make train-seed-789    Train with seed=789"
	@echo "  make train-all-seeds   Train all 5 seeds"
	@echo ""
	@echo "Validation:"
	@echo "  make analyze-convergence   5-model Procrustes analysis"
	@echo "  make viral-validation      15-virus κ sweeps"
	@echo "  make tree-validation       Phylogenetic tree κ estimation"
	@echo "  make validate-all          Run all validations"
	@echo ""
	@echo "Figures:"
	@echo "  make figures           Generate all publication figures"
	@echo ""
	@echo "Data:"
	@echo "  make fetch-data        Download public NCBI genomes"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean             Remove artifacts"

# =============================================================================
# Setup
# =============================================================================

install:
	@echo "📦 Installing dependencies..."
	pip install -r requirements.txt
	@echo "✅ Installation complete"

test:
	@echo "🧪 Running smoke tests..."
	python -c "from model import BiosphereCodec; import torch; \
		m = BiosphereCodec(100, 64, 1, 16); \
		l, _ = m(torch.randint(0, 99, (2, 32))); \
		print(f'✅ Model test passed (loss={l.item():.2f})')"

# =============================================================================
# Theory Verification
# =============================================================================

verify-lean:
	@echo "🔬 Verifying Lean proofs..."
	cd theory/lean && lake build
	@echo "✅ All proofs verified"

# =============================================================================
# Training
# =============================================================================

TRAIN_CMD = python model/training.py

train-seed-0:
	@echo "🧬 Training with seed=0..."
	$(TRAIN_CMD) --seed 0 --output_dir ./results/seed_0 --use_wandb False

train-seed-42:
	@echo "🧬 Training with seed=42..."
	$(TRAIN_CMD) --seed 42 --output_dir ./results/seed_42 --use_wandb False

train-seed-123:
	@echo "🧬 Training with seed=123..."
	$(TRAIN_CMD) --seed 123 --output_dir ./results/seed_123 --use_wandb False

train-seed-456:
	@echo "🧬 Training with seed=456..."
	$(TRAIN_CMD) --seed 456 --output_dir ./results/seed_456 --use_wandb False

train-seed-789:
	@echo "🧬 Training with seed=789..."
	$(TRAIN_CMD) --seed 789 --output_dir ./results/seed_789 --use_wandb False

train-all-seeds: train-seed-0 train-seed-42 train-seed-123 train-seed-456 train-seed-789
	@echo "✅ All 5 seeds trained"

# =============================================================================
# Validation
# =============================================================================

analyze-convergence:
	@echo "📊 Analyzing 5-model convergence..."
	python validation/genomic/scripts/5_model_convergence_analysis.py
	@echo "✅ Convergence analysis complete"

viral-validation:
	@echo "🦠 Running viral validation sweeps..."
	cd validation/viral/scripts && bash orchestration/run_all_validation_sweeps.sh
	@echo "✅ Viral validation complete"

tree-validation:
	@echo "🌳 Running phylogenetic tree validation..."
	python validation/phylogenetic/scripts/kappa_validation_pipeline.py \
		--files validation/phylogenetic/trees/*.nwk
	@echo "✅ Tree validation complete"

validate-all: verify-lean analyze-convergence tree-validation
	@echo "✅ All validations complete"

# =============================================================================
# Figures
# =============================================================================

figures:
	@echo "📈 Generating figures..."
	mkdir -p figures/outputs
	python figures/fig4_curvature_entropy.py
	@echo "✅ Figures saved to figures/outputs/"

# =============================================================================
# Data
# =============================================================================

fetch-data:
	@echo "📥 Fetching public NCBI genomes..."
	@echo "⚠️  This will download ~50GB of data"
	@read -p "Enter your email for NCBI: " email; \
	python data/scripts/fetch_from_manifest.py \
		--manifest data/manifests/public_refseq.tsv \
		--output ./raw_genomes \
		--email $$email
	@echo "✅ Download complete"

# =============================================================================
# Cleanup
# =============================================================================

clean:
	@echo "🧹 Cleaning..."
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
	rm -rf .lake theory/lean/.lake 2>/dev/null || true
	@echo "✅ Clean complete"
