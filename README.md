# Transcriptomics Pipeline (HISAT2 + StringTie + DESeq2)

This repository contains a reproducible pipeline for RNA-Seq Transcriptomic analysis using **HISAT2** for alignment and **StringTie** for transcript assembly and quantification, and **DESeq2** for differential gene expression analysis.

---

| Step | Description                                  | Tool                  |
| :--: | :------------------------------------------- | :-------------------- |
|  1️⃣ | **Quality Control**                          | FastQC                |
|  2️⃣ | **Trimming**                                 | Trim Galore           |
|  3️⃣ | **Alignment**                                | HISAT2                |
|  4️⃣ | **Transcript Assembly**                      | StringTie             |
|  5️⃣ | **Quantification & Count Matrix Generation** | StringTie, prepDE.py  |
|  6️⃣ | **Differential Expression Analysis**         | DESeq2                |
|  7️⃣ | **Visualization & Interpretation**           | R (ggplot2, pheatmap) |

---

## ⚙️ Environment Setup & Pipeline Execution

Before running the pipeline, create the required Conda environment using the provided YAML file:

```bash
conda env create -f env.yaml
```

Once the environment is successfully created, run the main script:
 using arguements: $1 is samples folder name,
                   $2 is Output folder,
                   $3 is Reference genome folder.

```bash
./Transcriptomics_pipeline.sh <reads_folder> <output_folder> <reference_genome>
```

---

## 🧫 Dataset Details

* **Organism:** *Homo sapiens* (GRCh38)
* **BioProjects:**

  * **Kidney:**

    * [PRJNA93899](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA93899) — 36 single-end samples
  * **Liver:**

    * [PRJNA750472](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA750472) — 12 single-end samples
    * [PRJNA429171](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA429171) — 6 paired-end samples

---

## 🔬 Biological Insights

* RNA-Seq analysis revealed distinct **tissue-specific expression** in Marburg virus–infected hosts.
* **Liver:** 1074 upregulated, 187 downregulated genes.
* **Kidney:** 972 upregulated, 363 downregulated genes.
* **Common** DEGs highlighted shared antiviral and immune responses.
* **Hub genes** (*ANXA1*, *ATF3*, *KLF4*, *TLR2*) showed key regulatory roles.
* *ATF3* and *ANXA1* interacted with viral proteins *VP30* and *VP35*, indicating immune modulation.
* **Lariciresinol** emerged as a potential inhibitor with good ADMET properties.
* Overall, the study identifies **key host targets and potential drug candidates** for Marburg virus infection.
