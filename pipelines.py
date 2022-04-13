import numpy as np
import sys, os
import pandas as pd
import scipy as sp
import scanpy as scp
from sklearn.decomposition import TruncatedSVD

class HiddenPrints:
    """
        Class to mute the print of a part of code.
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def quality_control_metrics(adata,mtGenesPositions):
    """
        Function that computes the basic quality control measures for the analysis and add them to the Ann data object.
    
        Input:
            adata: Annotated data in which compute the quality control measures.
            mtGenesPositions: Array of positions of the Mitochondrial Genes.
            
        Output:
            Annotated data of input with the following added information:
                - "#Counts": UMI count per cell added to obs.
                - "#Genes": Number of non-zero genes added to obs.
                - "mtFraction": Fraction of mitochondrial counts vs total added to obs.
                - "#Cells": Number of cells having non-zero counts of this gene. Added to var. 
    """
    adata.obs["#Counts"] = np.array(np.sum(adata.X,axis=1))[:,0]
    adata.obs["#Genes"] = np.array(np.sum(adata.X>0,axis=1))[:,0]
    adata.obs["mtFraction"] = np.array(np.sum(adata.X[:,mtGenesPositions],axis=1))[:,0]/adata.obs["#Counts"]
    adata.var["#Cells"] = np.array(np.sum(adata.X>0,axis=0))[0,:]
    
def scrublet_pipeline(adata,
                      batch_key=None,
                      normalize=True,
                      target_sum=None,
                      log_normalize=True,
                      feature_reduction_flavor="seurat",
                      n_retained_features=None,
                      exclude_genes=None,
                      key_exclude_list="Gene",
                      n_prin_components = 50,
                      knn_dist_metric="correlation",
                      key_added="Scrublet_score",
                      doublet_statistics_file=None,
                      verbose=False):
    """
        Function that customized the scrublet pipeline to compute it for different samples intependently and a personalised preprocessing.
        
        The function expects a raw count matrix.
    
        Input:
            adata:Annotated data in which compute the quality control measures.
            batch_key=None: Batch key over which perform the analysis independently. If None, the dataset is considered a unique sample.
            normalize=True: If normalizing the scores after simulating the doublets.
            target_sum=None: Target sum of the normalization.
            log_normalize=True: If log normalize the data.
            feature_reduction_flavor="seurat": Flavor of feature selection (see scanpy.pp.highly_variable_genes for possible options.)
            n_retained_features=None. Number of retained features, if none with "seurat", they are automatically selected.
            exclude_genes=None: Genes to be excluded even if selected as Higly Variable (cell cycle, sex genes...)
            key_exclude_list="Gene": Key by which the names of genes can be compared to the list of excluded_genes.
            n_prin_components = 50: Number of PCs to keep.
            knn_dist_metric="correlation": Metric of the KNN for computing during the execution of the scrublet algorithm.
            key_added="Scrublet_score": Key name added to the annotated data with the scrublet score.
            doublet_statistics_file=None: If different to None, a file is created with the scrublet score of the doublet statistics. This metric help to identify the threshold for imputation. 
            verbose=False: Verbosity of the scrublet algorithm.
        Output:
            Annotated data with the scrublet score added to .obs.
            Additionally, if doublet_statistics_file argument has been given, it creates a file with the simulated doublets score.
    """

    adata.obs[key_added] = 0
    
    if doublet_statistics_file != None:
        statistics = pd.DataFrame()

    if batch_key != None:
        for stage in adata.obs[batch_key].unique():
                
            #Extract stage
            b = adata[adata.obs[batch_key]==stage,:].copy()

            #Simulate the doublets
            bDoublets = scp.external.pp.scrublet_simulate_doublets(b,sim_doublet_ratio=2,random_seed=1)
                
            #Custom part of the pipeline
            if normalize:
                scp.pp.normalize_total(b,target_sum=target_sum)
                scp.pp.normalize_total(bDoublets,target_sum=target_sum)

            if log_normalize:
                scp.pp.log1p(b) #Logarithmize data
                scp.pp.log1p(bDoublets) #Logarithmize data
                    
            if feature_reduction_flavor != None:
                scp.pp.highly_variable_genes(b,flavor=feature_reduction_flavor) #Seurat flavour of obtaining HVGs
                    
                if exclude_genes != None:
                    b.var.loc[:,"higly_variable"] = [not ((g in exclude_genes) or not hv) for g,hv in b.var.loc[:,[key_exclude_list,"highly_variable"]].values]
                        
                bDoublets = bDoublets[:,b.var["highly_variable"]]
                b = b[:,b.var["highly_variable"]]

            #Compute the doublet factor
            if not verbose:
                with HiddenPrints():
                    scp.external.pp.scrublet(b,bDoublets,
                                                knn_dist_metric = knn_dist_metric,
                                                n_prin_comps = n_prin_components,
                                                log_transform=False)
            else:
                scp.external.pp.scrublet(b,bDoublets,
                                            knn_dist_metric = knn_dist_metric,
                                            n_prin_comps = n_prin_components,
                                            log_transform=False)
                
            if doublet_statistics_file != None:

                statistics_batch = pd.DataFrame()
                statistics_batch.loc[:,"doublet_score_sim"] = b.uns['scrublet']['doublet_scores_sim']
                statistics_batch.loc[:,batch_key] = stage
                statistics = statistics.append(statistics_batch).reset_index(drop=True)

            adata.obs.loc[b.obs.index,key_added] = b.obs.loc[:,"doublet_score"]
    else:
        #Extract stage
        b = adata.copy()

        #Simulate the doublets
        bDoublets = scp.external.pp.scrublet_simulate_doublets(b,sim_doublet_ratio=2,random_seed=1)
                
        #Custom part of the pipeline
        if normalize:
            scp.pp.normalize_total(b,target_sum=target_sum)
            scp.pp.normalize_total(bDoublets,target_sum=target_sum)

        if log_normalize:
            scp.pp.log1p(b) #Logarithmize data
            scp.pp.log1p(bDoublets) #Logarithmize data
                    
        if feature_reduction_flavor != None:
            scp.pp.highly_variable_genes(b,flavor=feature_reduction_flavor) #Seurat flavour of obtaining HVGs
                    
            if exclude_genes != None:
                b.var["higly_variable"] = [not ((g in exclude_genes) or not hv) for g,hv in b.var.loc[:,[key_exclude_list,"highly_variable"]].values]
                        
            bDoublets = bDoublets[:,b.var["highly_variable"]]
            b = b[:,b.var["highly_variable"]]

        #Compute the doublet factor
        if not verbose:
            with HiddenPrints():
                scp.external.pp.scrublet(b,bDoublets,
                                            knn_dist_metric = knn_dist_metric,
                                            n_prin_comps = n_prin_components,
                                            log_transform=False)
        else:
            scp.external.pp.scrublet(b,bDoublets,
                                        knn_dist_metric = knn_dist_metric,
                                        n_prin_comps = n_prin_components,
                                        log_transform=False)
                
        adata.obs.loc[b.obs.index,key_added] = b.obs.loc[:,"doublet_score"]
                
        if doublet_statistics_file != None:

            statistics_batch = pd.DataFrame()
            statistics_batch["doublet_score_sim"] = b.uns['scrublet']['doublet_scores_sim']
            statistics_batch[batch_key] = stage
            statistics = statistics.append(statistics_batch).reset_index(drop=True)

    if doublet_statistics_file != None:
        statistics.to_csv(doublet_statistics_file)

    return

def dimensionality_reduction_pipeline(adata,
                      add_key="X_pca_Stage",
                      batch_key="Stage",
                      feature_reduction_flavor="seurat",
                      n_retained_features=None,
                      exclude_genes=None,
                      key_exclude_list="Gene",
                      n_comps = 50):
    """
        Function that automatizes the dimensionality reduction pipeline over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            key_added_pca="X_pca_Stage": Key to be added to the pca components of the pipeline. To be added in obsm.
            batch_key="Stage": Batch key over which perform the analysis independently.
            feature_reduction_flavor="seurat": Flavor of feature selection (see scanpy.pp.highly_variable_genes for possible options.)
            n_retained_features=None. Number of retained features, if none with "seurat", they are automatically selected.
            exclude_genes=None: Genes to be excluded even if selected as Higly Variable (cell cycle, sex genes...)
            key_exclude_list="Gene": Key by which the names of genes can be compared to the list of excluded_genes.
            n_comps = 50: Number of PCs to keep.
        Output:
            Annotated data with the pca after the pipeline added to .obsm.
    """
    
    adata.obsm[add_key] = np.zeros([adata.shape[0],n_comps])

    for stage in adata.obs[batch_key].unique():
                        
        #Extract stage
        l = adata.obs[batch_key]==stage
        b = scp.AnnData(adata.X[l,:].copy())
        b.obs.index = adata.obs[l].index
                    
        #Compute Feature selection
        if feature_reduction_flavor != None:
            scp.pp.highly_variable_genes(b,flavor=feature_reduction_flavor) #Seurat flavour of obtaining HVGs
                    
            if exclude_genes != None:
                b.var.loc[:,"higly_variable"] = [not ((g in exclude_genes) or not hv) for g,hv in b.var.loc[:,[key_exclude_list,"highly_variable"]].values]
                    
        #Compute pca using TrucatedSVD
        X = TruncatedSVD(n_components=50,random_state=0).fit_transform(b.X[:,b.var.loc[:,"highly_variable"].values])
        
        adata.obsm[add_key][b.obs.index.values.astype(int),:] = X
        
    return

def harmony_integrate_batches(adata,
                      basis="X_pca_Stage",
                      add_key="X_pca_Harmony_Stage",
                      batch_key="Stage",
                      sample_key="Sample",
                      *args,**kwargs):
    """
        Function that automatizes the Harmony integration pipeline over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            basis="X_pca_Stage": Key of obsm defining the PCA components over which correct the basis.
            add_key="X_pca_Harmony_Stage": Key to be added to the pca components of the pipeline. To be added in obsm.
            batch_key="Stage": Batch key over which perform the analysis independently.
            sample_key="Sample": Key to specify the samples in the batch to be corrected.
            args: Arguments to be passed to harmony_integrate.
            kwargs: Keyword arguments to be passed to harmony_interagte.
        Output:
            Annotated data with the harmony corrected pca after the batch correction added to .obsm.
    """
    
    adata.obsm[add_key] = np.zeros_like(adata.obsm[basis])

    for stage in adata.obs[batch_key].unique():
                
        #Extract stage
        l = adata.obs[batch_key]==stage
        b = scp.adata(adata.obsm[basis][l,:])
        b.obs.index = adata.obs[l].index
        b.obs = adata.obs.loc[l,:]
        b.obs.loc[:,sample_key] = b.obs.loc[:,sample_key].astype(str)
                    
        #Compute harmony
        scp.external.pp.harmony_integrate(b,key=sample_key,basis=basis,*args,**kwargs)
        
        adata.obsm[add_key][b.obs.index.values.astype(int),:] = b.obsm["X_pca_harmony"]
        
    return

def neighbors_batches(adata,
                      use_rep="X_pca",
                      add_key="Neighbors_Stage",
                      batch_key="Stage",
                      *args,**kwargs):
    """
        Function that automatizes the KNN over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            use_rep="X_pca": Rep to use for the neighbors search. To be added in uns.
            add_key="Neighbors_Stage": Batch key over which perform the analysis independently. If None, the dataset is considered a unique sample.
            sample_key="Sample": Batch key over which perform the analysis independently.
            args: Arguments to be passed to harmony_integrate.
            kwargs: Keyword arguments to be passed to scp.pp.neighbors.
        Output:
            Annotated data with the neighbors by batches stored in uns.
    """

    adata.uns[add_key] = {}
    adata.uns[add_key]["connectivities"] = sp.sparse.lil_matrix((adata.shape[0],adata.shape[0]))
    adata.uns[add_key]["distances"] = sp.sparse.lil_matrix((adata.shape[0],adata.shape[0]))

    for stage in adata.obs[batch_key].unique():
        
        #Extract stage
        l = adata.obs[batch_key]==stage
        b = scp.AnnData(adata.obsm[use_rep][l,:].copy()) 
        b.obs.index = adata.obs[l].index
        
        scp.pp.neighbors(b,use_rep="X",*args,**kwargs)

        for i in b.uns["neighbors"].keys():
            adata.uns[add_key][i] = b.uns["neighbors"][i]

        l1 = b.obs.index.values.astype(int)[b.uns["neighbors"]["connectivities"].nonzero()[0]]
        l2 = b.obs.index.values.astype(int)[b.uns["neighbors"]["connectivities"].nonzero()[1]]
        v = np.array(b.uns["neighbors"]["connectivities"][b.uns["neighbors"]["connectivities"].nonzero()])[0,:]
        adata.uns[add_key]["connectivities"][l1,l2] = v

        l1 = b.obs.index.values.astype(int)[b.uns["neighbors"]["distances"].nonzero()[0]]
        l2 = b.obs.index.values.astype(int)[b.uns["neighbors"]["distances"].nonzero()[1]]
        v = np.array(b.uns["neighbors"]["distances"][b.uns["neighbors"]["distances"].nonzero()])[0,:]
        adata.uns[add_key]["distances"][l1,l2] = v
        
    adata.uns[add_key]["connectivities"] = adata.uns[add_key]["connectivities"].tocsr()
    adata.uns[add_key]["distances"] = adata.uns[add_key]["distances"].tocsr()
        
    return

def umap_batches(adata,
                      neighbors_key="Neighbors_Stage",
                      add_key="X_umap_Stage",
                      batch_key="Stage",
                      *args,**kwargs):
    """
        Function that automatizes the umap creation over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            basis="X_pca_Stage": Key of obsm defining the PCA components over which correct the basis.
            add_key="X_umap_Stage": Key to be added to the pca components of the pipeline. To be added in obsm.
            batch_key="Stage": Batch key over which perform the analysis independently.
            args: Arguments to be passed to harmony_integrate.
            kwargs: Keyword arguments to be passed to harmony_interagte.
        Output:
            Annotated data with the umap after the batch correction added to .obsm.
    """
    
    adata.obsm[add_key] = np.zeros([adata.shape[0],2])

    for stage in adata.obs[batch_key].unique():
                
        #Extract stage
        b = adata[adata.obs[batch_key]==stage,:].copy()
                    
        #Compute harmony
        scp.tl.umap(b,neighbors_key=neighbors_key,*args,**kwargs)
        
        adata.obsm[add_key][b.obs.index.values.astype(int),:] = b.obsm["X_umap"]
        
    return

def louvain_batches(adata,
                      neighbors_key="Neighbors_Stage",
                      add_key="Louvain_Stage",
                      batch_key="Stage",
                      resolution=.1,
                      *args,**kwargs):
    """
        Function that automatizes the louvain clustering over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            neighbors_key="Neighbors_Stage": Key of neighbors over which cluster the basis.
            add_key="Louvain_Stage": Key to be added with the clusters. To be added in obs.
            batch_key="Stage": Batch key over which perform the analysis independently.
            resolution=.1: Resolution of the clustering.
            args: Arguments to be passed to harmony_integrate.
            kwargs: Keyword arguments to be passed to harmony_interagte.
        Output:
            Annotated data with the harmony corrected pca after the batch correction added to .obsm.
    """
    
    adata.obsm[add_key] = np.zeros_like(adata.obsm[basis])

    for stage in adata.obs[batch_key].unique():
                
        #Extract stage
        b = adata[adata.obs[batch_key]==stage,:].copy()
        
        #Make clustering
        scp.tl.louvain(b,neighbors_key=neighbors_key,resolution=resolution,*args,**kwargs)
        
        #Add to the original data
        adata.obs.loc[b.obs.index.values,add_key] = b.obs.loc[:,"louvain"]
        
    return

def leiden_batches(adata,
                      neighbors_key="Neighbors_Stage",
                      add_key="Louvain_Stage",
                      batch_key="Stage",
                      resolution=.1,
                      *args,**kwargs):
    """
        Function that automatizes the leiden clustering over batches.
        
        Input:
            adata:Annotated data in which compute the quality control measures.
            neighbors_key="Neighbors_Stage": Key of neighbors over which cluster the basis.
            add_key="Louvain_Stage": Key to be added with the clusters. To be added in obs.
            batch_key="Stage": Batch key over which perform the analysis independently.
            resolution=.1: Resolution of the clustering.
            args: Arguments to be passed to harmony_integrate.
            kwargs: Keyword arguments to be passed to harmony_interagte.
        Output:
            Annotated data with the harmony corrected pca after the batch correction added to .obsm.
    """
    
    adata.obsm[add_key] = np.zeros_like(adata.obsm[basis])

    for stage in adata.obs[batch_key].unique():
                
        #Extract stage
        b = adata[adata.obs[batch_key]==stage,:].copy()
        
        #Make clustering
        scp.tl.leiden(b,neighbors_key=neighbors_key,resolution=resolution,*args,**kwargs)
        
        #Add to the original data
        adata.obs.loc[b.obs.index.values,add_key] = b.obs.loc[:,"louvain"]
        
    return