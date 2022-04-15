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

In this file, we perform clustering and annotation of the different stages of the dataset.

### Clustering

Following the same reasoning as previously stated, we perfrom the clustering and annotation over each time point independently.

We cluster the different stages using a the following algorithms:

1. Louvain

In order to check if the clusters are biologically relevant, we check the distribution of doublet expression in those clusters. We can see that several clusters contain a specially high level of doublet expression in comparison with the average. Marking those clusters that clusters that are 2 standard deviations over the average doublet score per cluster, we see that the majority of those clusters are small populations of outliers or belong to regions between clusters, as can be clearly seen in stages "E8.25" and "E8.5". The presence of these clusters with atypical overexpression of doublet scores may indicate that we did not remove all possible doublets from the dataset.

### Annotation

The manual annotation is performed over each stage independently.

For the manual annotation we use the following information:

1. A list of Differentially Expressed genes of each cluster vs. the rest.
2. Plots of gene expression over the UMAP for characteristic genes (list of genes to be found in parameters.py).

The dataset was previously annotated by the authors of the dataset. In order to check the potential discrepancies between annotations:

 1. We make comparative plots of the annotations (each column adds to one).
 2. Plots of their clusters and annotations over our representation.

 ## 3 Clustering and Annotation.ipynb

 ### Dimensionality reduction of the full dataset

Until now, the analysis has been performed at each stage independently. Now we are interested in making a full integration of the data to see the transcriptomic trajectories along time.

For that, we perform a dimensionality reduction using the following pipeline:

 1. Selection of Highly Varying Genes using `seurat` flavor with the selection being the combination of the HVGs of each stage. We remove genes related with sex and cell cycle from this list.
 2. PCA, in the flavor of Truncated SD Decomposition.
 
Notice that for this embedding we do not make any kind of batch correction as it is not clear what should be the correction between stages and it does not make much sense to make a correction between samples of each batch while mantaining the rest fixed.
 
This returns us a latent space of the hole dataset.

### Paga by stages and samples
We expect that the connections between cells are more strongly connected in a causal order, being stages at close points having cells more connected than others.

We can test this making a PAGA graph. For that we perform the following process:

 1. Make a KNN neighbour graph over the reduced space mentioned above.
 2. Make PAGA graph by Stages and by Samples.


We can see that the stages are connected in a time order with a few exceptions:

 1. Mixed gastrulation is considered an outlier. This makes sense as it is a combination of cell types. We can conclude that this subset of the data is not informative for the analysis that we are performing.
 2. Stages E6.5 to E7.25 are very connected between them. There can be two reasons behind this exploration. the first is that transcriptomically the development between these stages is slower, making them very close to each other. The second reason is that the different embryos composing the different samples of each embryo are not annotated in the correcponding stage with high precission, giving rise to Stage annotation missalignements that lead to the high clustering. This second reason seems to be backed when the PAGA analysis is performed over the sample graph.
 3. A similar thing happen between stages E7.75 to E8.25.

 We finally, we in a UMAP the integrated dataset. The complexity and density of the data does not allow to see much.

 ### Transcriptomical connection between stages

We connect the clusters between stages in fashion very similar to that proposed in [this paper](https://www.nature.com/articles/s41588-022-01018-x).

The connections are made in to manners:

 1. Forward: Gives information of the clusters in a earlier stage to which future clusters provide in the future.
 2. Backward: Gives information of where the clusters in a later stage comes from.