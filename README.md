# DeepChrome & AttentiveChrome — PyTorch Replication (Gene expression prediction from Histone modifications) 

A PyTorch replication of [DeepChrome](https://academic.oup.com/bioinformatics/article/32/17/i639/2450757?searchresult=1#394406142) (Singh et al., 2016) and its attention-based extension [AttentiveChrome](https://arxiv.org/abs/1708.00336) (Singh et al., 2017), trained and evaluated on cell type E047 (Primary T CD8+ naive cells) obtainedd from the NIH Roadmap Epigenomics Mapping Consortium (REMC).

---

## Project Overview

DeepChrome predicts gene expression (high/low) from histone modification signals using a CNN. AttentiveChrome extends this with a hierarchical LSTM + attention mechanism, adding biological interpretability on top of competitive predictive performance.

This project replicates both models from scratch in modern PyTorch, starting from the original Lua/Torch7 source, and evaluates them on E047.

---

## Repository Structure

```
deepchrome/
│
├── Preprocessing
│   ├── extract_tss.py               # Extract TSS coordinates from GENCODE v19 GTF
│   ├── create_bins.py               # Divide ±5kb TSS windows into 100bp bins
│   ├── convert_to_bam.py            # Convert REMC tagAlign files to BAM format
│   └── count_reads.py               # Count histone mark reads per bin (bedtools multicov)
│
├── Training & Evaluation
│   ├── deepchrome_1_preprocessing.ipynb     # Dataset assembly and feature engineering
│   ├── deepchrome_2_authers_model.ipynb     # DeepChrome model, training, Optuna tuning
│   └── AttentiveChrome_Replication.ipynb    # AttentiveChrome model, training, attention visualization
```

---

## Dataset

- **Cell type**: E047 — Primary T CD8+ naive cells (REMC)
- **Histone marks**: H3K4me3, H3K4me1, H3K36me3, H3K9me3, H3K27me3
- **Genes**: 19,300 after filtering (GENCODE v19, chromosomes 1–22 + X)
- **Features**: 500 per gene — 100 bins × 5 histone marks, read counts aggregated per 100bp bin in a ±5kb TSS window
- **Labels**: Binary — high (1) / low (0) expression based on RPKM median split
- **Split**: Train / Val / Test — 6,601 / 6,601 / 6,098 (`np.random.seed(1)`)

---

## Models

### DeepChrome

A CNN that takes histone modification signals as a (5 marks × 100 bins) feature map.

```
Conv1d(5 → 50, kernel=10) → ReLU → MaxPool1d(5)
→ Dropout → Linear(900 → 625) → ReLU
→ Linear(625 → 125) → ReLU → Linear(125 → 2)
```

- **Loss**: NLLLoss (with log_softmax)
- **Optimizer**: SGD (lr=1e-3, momentum=0)
- **Batch size**: 1

### AttentiveChrome

A hierarchical LSTM + attention model that encodes each histone mark's bin sequence independently, then attends across marks.

**Three components:**
1. **BinEncoder** — Bidirectional LSTM encoding bin sequences per mark
2. **Attention (α)** — Soft attention over bins for each mark (bin-level)
3. **Attention (β)** — Soft attention over marks (mark-level)

```
Input (N, 100, 5) → 5 × BinEncoder (BiLSTM) → 5 × Attention(α)
→ Mark representations → BinEncoder (BiLSTM) → Attention(β)
→ Linear(1) → Sigmoid
```

- **Loss**: BCELoss
- **Optimizer**: Adam (lr=1e-4)
- **Batch size**: 32

---

## Results

| Configuration | Test AUROC |
|---|---|
| DeepChrome — paper hyperparameters, 100 epochs | 0.9061 |
| DeepChrome — early stopping (patience=10) | 0.9189 |
| DeepChrome — Optuna tuning + early stopping | 0.9201 |
| AttentiveChrome — early stopping (patience=10) | 0.9198 |
| **Paper reported average (56 cell types)** | **0.8000** |

Both models substantially exceed the paper's reported average AUROC of 0.80. AttentiveChrome matches DeepChrome Optuna performance without any hyperparameter tuning, while also providing biological interpretability through its attention weights.

Early stopping was the single most impactful regularization step for both models. Optuna tuning on DeepChrome confirmed the original authors' hyperparameters were already near-optimal, with only marginal gain (0.9189 → 0.9201).

---

## Attention Visualization — E047

Attention weights were averaged across all test set genes to produce a single representative view of what the model learned to focus on.

![Attention Visualization](attentivechrome_attention_e047.png)

### Mark-level Attention (β)

| Mark | β Weight | Role |
|---|---|---|
| H3K4me3 | 0.345 | Active promoter mark |
| H3K27me3 | 0.221 | Polycomb repression |
| H3K4me1 | 0.177 | Enhancer mark |
| H3K9me3 | 0.137 | Constitutive heterochromatin |
| H3K36me3 | 0.120 | Transcribed gene body |

**H3K4me3 dominance** is expected — it is the canonical active promoter mark, tightly focused at the TSS, and the strongest predictor of gene expression. Its high β weight reflects its direct relevance to transcriptional activation in CD8+ naive T cells.

**H3K27me3 as second highest (β=0.221)** is biologically meaningful. In naive CD8+ T cells, bivalent chromatin domains — where H3K4me3 and H3K27me3 co-occur — are well documented at genes that are silenced but poised for rapid activation upon antigen encounter. The model's joint weighting of both marks likely captures this bivalency, where the co-occurrence pattern is more informative than either mark alone.

**H3K4me1's moderate weight (β=0.177)** suggests the model is also capturing distal regulatory activity via enhancers, which are marked by H3K4me1 rather than H3K4me3.

### Bin-level Attention (α)

**H3K4me3** shows the sharpest attention peak at the TSS (bin 50), consistent with its well-characterized tight enrichment at active promoters.

**H3K36me3** shows the most distinctive spatial pattern: low attention near the TSS and elevated attention toward the distal ends of the window (bins 0–5 and ~95–100). This reflects its known biology — H3K36me3 is deposited by SETD2 during transcriptional elongation and accumulates in gene bodies downstream of the TSS, not at the promoter itself. The model learned this spatial distribution without explicit positional supervision.

**H3K27me3 and H3K9me3** show relatively uniform attention across bins, consistent with their role as broad repressive domains not tightly anchored to the TSS.

---

## Usage

### 1. Preprocessing (Linux VM)

Run the preprocessing scripts in order to build the feature matrix from raw REMC data:

```bash
python extract_tss.py        # Extract TSS from GENCODE v19 GTF
python create_bins.py        # Generate 100bp bins per TSS window
python convert_to_bam.py     # Convert tagAlign files to BAM
python count_reads.py        # Count reads per bin per mark
```

The final dataset (`E047_dataset.csv`, 19,300 × 501) is assembled in `deepchrome_1_preprocessing.ipynb`.

### 2. DeepChrome Training & Evaluation

Open `deepchrome_2_authers_model.ipynb` on Kaggle (GPU recommended). The notebook covers:
- Model definition (faithful to original Lua architecture)
- Training with early stopping
- Optuna hyperparameter tuning (30 trials)
- Test set evaluation

### 3. AttentiveChrome Training & Evaluation

Open `AttentiveChrome_Replication.ipynb` on Kaggle (GPU recommended). The notebook covers:
- Model definition (BinEncoder, Attention, AttentiveChrome)
- Training with early stopping
- Test set evaluation
- Attention weight extraction and visualization

---

## Requirements

```
torch
numpy
pandas
scikit-learn
optuna
matplotlib
bedtools (for preprocessing)
samtools (for preprocessing)
```

---

## References

- Singh, R., Lanchantin, J., Robins, G., & Qi, Y. (2016). [DeepChrome: Deep-learning for predicting gene expression from histone modifications](https://arxiv.org/abs/1607.02078). *Bioinformatics*.
- Singh, R., Lanchantin, J., Sekhon, A., & Qi, Y. (2017). [AttentiveChrome: Attend and Predict: Understanding Gene Expression with Deep Neural Networks](https://arxiv.org/abs/1708.00336). *NeurIPS*.
- NIH Roadmap Epigenomics Mapping Consortium. https://www.roadmapepigenomics.org/
- GENCODE v19. https://www.gencodegenes.org/human/release_19.html
