# Analysis of Pijuan-Sala dataset

This dataset comes from the paper of [Pijuan-Sala et al. (2019)](https://www.nature.com/articles/s41586-019-0933-9).

For this analysis, we depart from the count matrix obtained after the alignement of the reads and do a systematic analysis.

The structure of the analysis is divided in different scrips that are described in the following.

## 0 Construct AnnData.ipynb

From the count matrix provided by the original authors, we create the Annotated Data object for the rest of the pipeline.

## 1 QC.ipynb

In here, we perform a Quality Control of the data. 

The quality control is based on the following metrics:

1. Total count of reads per cell
2. Total count of non-zero genes per cell
3. Mitochondrial fraction

The dataset as already provided by the authors has removed the tails of the total counts. The only group of cells that have not been removed are stripped cells, as detected by the bimodality at very low mitochondrial fractions. We set a threshold for imputing these cells aas sttripped cells.

Before this control, we account for doublet cells by running a Scrublet algorithm for the detectiuon of doublets. We perform this algorithm with the following characteristics:

- Each sample is analized independently as doublets are produced within each mouse sample with their own batch effects and disgregation problems.
- The doublet score is computed generating doublets and considering the neighbors in a reduced space constructed as:
    - Highly Variable Genes (HVGs) as selected automatically by "seurat" algorithm in `scanpy.pp.highly_varying_genes`.
    - Removed HVGs related with cell cycle and sex.
    - Retained 50 Principal Components (PCs) `scanpy.pp.pca`.
    - KNN matrix generated with the 50 PCs and "correlation" metric as implemented in `scanpy.pp.neighbors`.

We choose the doublet score to impute as cells base on a unique threshold for all doublet bimodality diagrams. (To be done: Performing an alternative control performing the doublet score trhough an alternative method does not change scores obtained.)

After this procedure, we compare the number of cells imputed by our procedure and the cells removed during the analysis of the original paper. We notice a big discrepancy of imputed cells, being our procedure far more permissive.

Finally, we removed the imputed cells from the dataset and save the new data.

## 2 Dimensionality Reduction.ipynb

