# Analysis of Pijuan-Sala dataset (v0.1)

This dataset comes from the paper of [Pijuan-Sala et al. (2019)](https://www.nature.com/articles/s41586-019-0933-9).

For this analysis, we depart from the count matrix obtained after the alignement of the reads and do a systematic analysis.

The structure of the analysis is divided in different scrips that are described in the following.

## 0 Construct AnnData.ipynb

From the count matrix provided by the original authors, we create the Annotated Data object for the rest of the pipeline.

## 1 Quality Control and Dimensionality Reduction.ipynb

### Quality Control

In here, we perform a Quality Control of the data. 

The quality control is based on the following metrics:

1. Total count of reads per cell
2. Total count of non-zero genes per cell
3. Mitochondrial fraction

The dataset as already provided by the authors has removed the tails of the total counts. The only group of cells that have not been removed are stripped cells, as detected by the bimodality at very low mitochondrial fractions. We set a threshold for imputing these cells as stripped cells.

The scatter plot shows an anomalous tendency formed by a second mode on the distribution of the data. This second modality is comming from hematopoietic cells. These cells have a lower varibility of genes and read counts as the nucleus is shutted down, which explains the atypical count number and the origin of the bimodality.

Before this control, we account for doublet cells by running a Scrublet algorithm for the detectiuon of doublets. We perform this algorithm with the following characteristics:

- Each sample is analized independently as doublets are produced within each mouse sample with their own batch effects and disgregation problems.
- The doublet score is computed generating doublets and considering the neighbors in a reduced space constructed as:
    - Highly Variable Genes (HVGs) as selected automatically by "seurat" algorithm in `scanpy.pp.highly_varying_genes`.
    - Removed HVGs related with cell cycle and sex.
    - Retained 50 Principal Components (PCs) `scanpy.pp.pca`.
    - KNN matrix generated with the 50 PCs and "correlation" metric as implemented in `scanpy.pp.neighbors`.

We choose the doublet score to impute as cells base on a unique threshold for all doublet bimodality diagrams. Up to this point of the analysis we do not remove any cells from the dataset.
### Dimensionality reduction

In this stage, we make a dimensionality reduction of the RNA space that encodes the information in a relevant space. **This analysis is performed for each temporal Stage of development independently as different times during development may be captured more effectively with different effective spaces.**

For obtaining each reduced space we follow the following pipeline:

1. Normalize the total count to the mean number of counts among all cells of the dataset.
2. Log normalize the counts. 
3. Select Highly Varying Genes using Seurat algorithm as implemented by `scanpy.pp.highly_varying_genes`. We remove from the highly varying genes those genes related to growth and sex.
4. Perform PCA analysis over the HVGs divemnsions and retain the 50 most variable PCs as our reduced space. We use the implementation of `scanpy.pp.pca`.

For visualization of the data, we further perform the following steps:

5. KNN matrix generated with the 50 PCs and "correlation" metric as implemented in `scanpy.pp.neighbors`.
6. UMAP in 2D representations as implemented in `scanpy.tl.umap`.

Visualization of the reduced space by stages shows batch effects distortions that should be corrected. We correct the batch effects with the follwing algorithms:

1. Harmony provided in `scanpy.external.pp.harmony_integrate`.

The PCS corrected show a proper correction of the batch effects between samples for all stages.

### Imputed cells

After correcting for batch effects, we visualize the location of stripped and doublets cells that where imputed during the Quality Control. Most of the imputed cells show grouped in clusters. This is consistent with our expectation of stripped and doublet cells to group as "cell types" of artificial nature. 

For comparison, we visualize the cells that we impute in our analysis and those from the original paper from Marioni finding a great discrepancy in the removal of cells. Putting a more restrictive threshold to our analysis recover a closer cleaning of the original analysis showing that our filtering is in general more permissive.

Having inspected the reasoable outcome of the imputed cells, we proceed to remove those and recompute the reduced space of the cleaned dataset from steps 1-6 described above and save the data.

## 2 Clustering and Annotation.ipynb

