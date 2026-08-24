# q2-eco-processes: Official QIIME 2 Plugin for Microbial Ecological Assembly Processes

[![Bioconda](https://img.shields.io/badge/bioconda-q2--eco--processes-blue.svg)](https://bioconda.github.io/recipes/q2-eco-processes/README.html)
[![QIIME 2 Library](https://img.shields.io/badge/QIIME%202%20Library-eco--processes-purple.svg)](https://library.qiime2.org)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

> **Quantify Selection, Dispersal Limitation, and Drift in Microbial Communities**  
> Developed by **LaBiOmics / UMC (Bioinformatics & Microbiome Unit)**

---

## 🔬 Overview

`q2-eco-processes` is an official QIIME 2 plugin that quantifies the relative influence of the **5 Stegen Ecological Assembly Processes** (Stegen et al. 2013, 2015) and **Ning's Normalized Stochasticity Ratio (NST)** (Ning et al. 2019) using phylogenetic null models ($\beta\text{NTI}$ and $RC_{\text{bray}}$):

1. **`Homogeneous_Selection`** ($\beta\text{NTI} < -2$) — Strong environmental filtering selects similar taxa.
2. **`Variable_Selection`** ($\beta\text{NTI} > +2$) — Divergent environmental conditions select different specialized taxa.
3. **`Dispersal_Limitation`** ($|\beta\text{NTI}| \le 2$ and $RC_{\text{bray}} > +0.95$) — Physical or physiological barriers prevent dispersal.
4. **`Homogenizing_Dispersal`** ($|\beta\text{NTI}| \le 2$ and $RC_{\text{bray}} < -0.95$) — High dispersal rates swamp community differences.
5. **`Undominated_Drift`** ($|\beta\text{NTI}| \le 2$ and $|RC_{\text{bray}}| \le 0.95$) — Stochastic ecological drift, weak selection, and birth-death fluctuations.

---

## 📋 Best Practices & Methodological Recommendations

To ensure maximum statistical validity and avoid sampling or tree-rooting biases when using `q2-eco-processes`:

1. **Rarefy Feature Tables First:**  
   Always run `qiime feature-table rarefy` to normalize sequencing depth across samples before calculating $RC_{\text{bray}}$ distances.
2. **Filter Ultra-Rare ASVs:**  
   Use `--p-min-frequency 10` to remove low-frequency PCR/sequencing artifact ASVs before null model tree simulations.
3. **High-Quality Rooted Phylogeny:**  
   Generate your `Phylogeny[Rooted]` input using MAFFT + FastTree (`qiime phylogeny align-to-tree-mafft-fasttree`) to ensure reliable $\beta\text{MNTD}$ distances.
4. **Test Phylogenetic Signal:**  
   Use `qiime eco-processes check-phylogenetic-signal` to confirm that phylogenetic distance correlates with environmental niche conservatism.

---

## 💻 Installation

### Option 1: Via Conda / Bioconda (Recommended)

```bash
conda install -c bioconda -c conda-forge q2-eco-processes
```

### Option 2: Via PyPI

```bash
pip install q2-eco-processes
```

---

## 🚀 QIIME 2 CLI Usage Example

### 1. Main Assembly Quantification & Dashboard
```bash
# 1. Rarefy feature table to uniform sampling depth
qiime feature-table rarefy \
  --i-table table-filtered.qza \
  --p-sampling-depth 10000 \
  --o-rarefied-table table-rarefied.qza

# 2. Quantify ecological assembly processes across experimental groups (e.g. Control vs Treated)
qiime eco-processes calculate-processes \
  --i-table table-rarefied.qza \
  --i-tree rooted_tree.qza \
  --m-metadata-file metadata.tsv \
  --m-column Experimental_Group \
  --p-min-frequency 10 \
  --p-n-jobs 4 \
  --o-process-table process_table.qza

# 3. Summarize assembly profile in interactive dashboard (.qzv)
qiime eco-processes summarize-processes \
  --i-table process_table.qza \
  --o-visualization process_summary.qzv
```

### 2. Export Pairwise Distance Matrices ($\beta\text{NTI}$ and $RC_{\text{bray}}$) for PCoA / PERMANOVA
```bash
qiime eco-processes calculate-matrices \
  --i-table table-rarefied.qza \
  --i-tree rooted_tree.qza \
  --p-n-jobs 4 \
  --o-beta-nti beta_nti_matrix.qza \
  --o-rc-bray rc_bray_matrix.qza
```

### 3. Test Phylogenetic Signal (Mantel Test)
```bash
qiime eco-processes check-phylogenetic-signal \
  --i-table table-rarefied.qza \
  --i-tree rooted_tree.qza \
  --m-metadata-file metadata.tsv \
  --m-column pH \
  --o-visualization phylo_signal.qzv
```

---

## 📜 Citations

1. Stegen JC et al. (2013) Quantifying community assembly processes and identifying features that impose them. *The ISME Journal* 7(11):2069-2079.
2. Stegen JC et al. (2015) Estimating and mapping ecological processes influencing microbial community assembly. *Frontiers in Microbiology* 6:370.
3. Ning D et al. (2019) A general framework for quantitatively assessing ecological stochasticity. *PNAS* 116(34):16892-16898.
