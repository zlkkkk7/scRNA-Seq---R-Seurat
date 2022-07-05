# scRNA-Seq---R-Seurat
# Cluster genes

1. Load libraries 
    1. install packages: seurat, tidyverse [package - install - search (choose **CRAN**) - install]
    2. load library 
        
        ```r
        #libraries 
        library (Seurat)
        library (SeuratDisk)
        ```
        

1. setup the **seurat object**
    1. read 
        
        ```r
        #Load the dataset
        nsclc <- Read10X_h5(filename = "C:/Users/vivia/Downloads/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
        #check the data set: because there are multiplen modalities 
        str(nsclc) 
        #select one modalities (this is a sparse matrix)
        cts <- pbmc.data$Gene Expression
        ```
        
        - read different single cell data formats [https://www.youtube.com/watch?v=3xcTpqQzUwQ](https://www.youtube.com/watch?v=3xcTpqQzUwQ)
            
            ```r
            #R Data Format (.rds)
            readRDS('xxx.rds')
            
            #10x CellRanger (.hdf5)
            Read10X_h5(filename = "xxx.h5", 
            use.names=TRUE,
            unique.features = TRUE)
            
            #AnnData Object (.h5ad)
            
            #Loom (.loom)
            
            #test based Market Exchange Format (.mtx)
            ```
            
        - notes:
            - sparse matrix: matrix that is mostly composed on 0
    2. load to seurat object 
        
        ```r
        #function that used to create a seurat object
        CreateSeuratObject (counts=, project= "", min.cells= , min.features= )
        #create 
        nsclc.seurat.obj <- CreateSeuratObject (counts=cts , project= "NSCLC", min.cells=3 , min.features=200 )
        #check the seurat object 
        str(nsclc.seurat.obj)
        ```
        
        - what does each parameter in this function mean?
            - counts: the dataset that we want to load to seurat object
            - project: name the project
            - min.cells: save features that present in # cells
            - min.features: save genes that have # features
        - **what is a seurat object and what does each slot mean**
        
2. quality control 
    1. Create matrix
        1. matrix 1 [low: bad quality of cells; high: doublets, multiplets]
            1. number of genes in a cell 
            2. number of molecules in a cell
            
            already present in the “meta.data” slot
            
            ```r
            View(nsclc.seurat.obj@meta.data)
            ```
            
        2. matrix 2: mitochondrial gene %
            
            need to calculate 
            
            ```r
            #function used
            	#PercentageFeatureSet (nsclc.seurat.obj, pattern="^MT-")
            #add a column "percent.mt" in the seurat oobject
            nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet (nsclc.seurat.obj, pattern="^MT-")
            ```
            
        
        visualize the all the matrix
        
        ```r
        #violin plots
        VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
        #feature scatterplot and a regression line
        FeatureScatter(nsclc.seurat.obj, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") + geom_smooth(method = 'lm') 
        ```
        
    2. Filter 
        
        ```r
        #function subset()
        nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset= nFeature_RNA>200, & nFeatureRNA<2500 & percent.mt<5) 
        ```
        

1. Normalize data: to compare sequencing data between cells 
    
    ```r
    nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
    ```
    

1. Identify variable features
    
    ```r
    #Calculate 
    nsclc.seurat.obj <- FindVariableFeatures (nsclc.seurat.obj, selection.method="vst", nfeatures=2000)
    #identify 
    top10 <- head(VariableFeatures(nscls.seurat.obj),10)
    
    #plot 
    plot1 <- VariableFeaturePlot (nsclc.seurat.obj)
    LabelPoints (plot=plot1, points=top10)
    ####Warning messages:
    #1: Transformation introduced infinite values in continuous x-axis 
    #2: Removed 2954 rows containing missing values (geom_point).Transformation introduced infinite values in continuous x-axis 
    ```
    

1. Scaling: clear unwanted noise during clustering 
    
    ```r
    #pick all features
    all.genes <- rownames (nsclc.seurat.obj)
    #scale them 
    nsclc.seurat.obj <- ScaleData (nsclc.seurat.obj, features=all.genes)
    ```
    
    slots:
    
    assay slot
    
    RNA assay
    
    counts: raw
    
    data: log normalized
    
    scale.data: scaled data 
    
2. linear dimension reduction: source of heterogeneity in the sample
    
    ```r
    #reduction 
    nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, feature=VariableFeatures(object=nsclc.seurat.obj))
    
    #print the top 5 PC, 5 feature
    print(nsclc.seurat.obj [["pca"]], dim=1:5, nfeatures=5)
    #visualize in heatmap (each PC scores and features)
    DimHeatmap(nsclc.seurat.obj, dim=1, cells=500, balanced=TRUE)
    ```
    
    - dim: # of PC
    - nfeatures: # of features
- determine dimensionality
    - many PC, but pick ones that have higher variance (the more the better)
    
    ```r
    ElbowPlot(nsclc.seurat.obj)
    #after viewing the plot, 15 PC is reasonable 
    ```
    

8. clustering: cluster cells that have similar features profiles
    
    ```r
    #picked the first 15 PC
    nsclc.seurat.obj <- FindNeightbors(nsclc.seurat.obj, dim=1:15) 
    
    #place cells to each cluster, and try out different resolutions 
    nsclc.seurat.obj <- FindCluster(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
    
    #visualize
    View(nsclc.seurat.obj@meta.data) 
    DimPlot(nsclc.seurat.obj, group.by ="RNA_snn_res.0.1", label=TRUE)
    
    #assign one resolution 
    Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
    ```
    
- granualarity/resolution of the cluster: lower number, fewer clusters → try out different resolution to see which clustering works the best
- further view data into clusters in low dimension space
