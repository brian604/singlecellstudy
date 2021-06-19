### Singlecellstudy
Project (Soongsil University Graduate Study) and TIL(Today I Learned) on single cell

## COVID project 

- Author: Jiwon Kim 
- For Soongsil University Single Cell sequencing analysis graduate course (under Prof. Junil Kim) 
- Paper/ Data Reference: <br>
[Single-cell landscampe of bronchoalveolar immune cells in patients with COVID-19](https://www.nature.com/articles/s41591-020-0901-9)
- Data that is ready to be analyzed is in data/ in h5 format
	- TCR-seq is out of scope in this analysis report.
	- A normal subject (GSM3660650) has not been included in the analysis.
	 Meta data containing patient and barcode information is in .txt file, named all.cell.annotation.meta.txt or meta.txt 
- Metadata containing patient (stratification into healthy/moderate/severe) and barcode are in the data/ folder. 
- Folder: image directory is broken down into the followings: 
	- pre_qc : Pre-QCed samples 
	- post_qc: Post-QCed samples 

#### Introduction of scRNA sequencing analysis
______________

### Introduction 

<p> <t> COVID-19 is a global pandemic known to be originated from mainland China. 
    SARS-CoV2 virus is known to infect respiratory tract, first. Infected individuals have exhibited a large range of symptoms. Those symptoms are known to point to differential immune response (Paces <em>et al.</em> 2020). 
    However, a high-resolution respiratory immune landscape is largely unknown. Since brochoalveolar alvage fluid (BAL) mirrors local immune landscape, they have attempted to do single-cell sequencing on it to three subject populations: healthy control, moderately illed, and severely- illed. With the carefully laid-out definitions of severeity, they have sequenced BAL of a total of 13 patients including 3 healthy control subjects. 

### Method 
#### QC & Integration
________________
- QC was done accordingly to the paper 
    - Minimum RNA Features = 200
    - Maximum RNA Features = 6000
    - Required counts of RNA = 1000
    - Maximum mito. cut-off = 10
- Integration was also done accordingly to the paper (first 50 dimensions)
<p>
    
```R
# Of which have been QC(filtered), 
all <- c(healthy.df.filtered, moderate.df.filtered, severe.df.filtered)
nCoV <- FindIntegrationAnchors(object.list = all, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))
```

#### Clustering 
________________
- Clustering was done accordingly to the paper 
    - Normalization using 'LogNormalize' method 
    - 'vst' method to identify top 2000 variable genes 
    - Scaling was done with variables 'nCount_RNA' and 'percent'mito'.


### Results: 
_______
<br>

##### Observation of Pre-QC
______________
<p> After downloading the files (h5 format), I have looked at the raw distribution of (1) number of RNA features (2) number of RNA (reads) counts (3) percentage of mitochondria across samples. </p>
I have noticed the followings: <p> 
Normal people (corresponding to the C51, C52, and C100) all have distinct bimodal distribution of RNA feature numbers. On the other hands, COVID-19 infected people have variable distributions. Although it would be nice to do statistical test on the distribution shape (manifold?), I have not done the analysis. <p>
For three moderately infected patients (C141, C142, C144), I have noted that highly diverse distributions of RNA features; from bimodal to long-tailed distribution really close to 0. 
The rest of samples (N = 6) are severely illed patients. Note the general skewed distribution close to 0. Also, RNAs (in terms of both frequencies and kinds) are very diverse, observing from a lot of points concentrated at the bottom for severely illed patients.  </p>

**Figure 1.1: Pre-QC of representative normal subject (C51)**
![Figure 1.1: Pre-QC of representative normal subject (C51)](images/pre_qc/C51_qc.png)
**Figure 1.2: Pre-QC of representative moderatly illed subject (C142)**
![Figure 1.2: Pre-QC of representative moderately illed subject (C142)](images/pre_qc/C142_qc.png)
**Figure 1.3: Pre-QC of representative severely illed subject (C145)**
![Figure 1.3: Pre-QC of representative severely illed subject (C145)](images/pre_qc/C145_qc.png) 



### Observation of Integration and Post-QC 
<p> After post-QC and integration, a landscape of RNA features and read counts across subjects has changed greatly. 
    Distinct bi-modal distribution of RNA feature numbers has disappeared from the normal subjects data. Except for one subject (i.e. C100), fairly uniform distributions for both RNA feature numbers and frequencies were observed. On the other hands, bi-modal distribution was seen on moderately illed patients for both RNA feature numbers and frequencies (except for one subject). Increased diversities had been observed in severely illed patients. </p>
    <br> 

**Figure 2: Post-QC of all subjects**
![Figure 2: Post-QC of all subjects](images/post_qc/qc.png)

### Identifying distinct distribution patterns upon patient group (i.e. severeity) and COVID-19 presence 
<br> 
<p>
    Using the metadata to differentiate and identify clusters has failed due to several reasons: (1) barcode number mismatch (2) technical difficulties. Jumping into conclusion, resolving subject information into Seurat object has been successful. 
</br> </p>
<p>
    Hence, I have assessed the UMAP projection of all cells to look at the distinct clusters of which specfic subgroups (i.e. COVID-19 and non-COVID19). 
Figure 3 and Figure 4 have shown the distinct patterns do patients have. To delve more to see what clusters are specific to the subgroup, a proportion of clusters upon disease presence (Figure 5) has identified which clusters deserve more attention. Clusters have been selected upon (1) relative contribution of clusters that form different disease groups, and (2) modest convolutedness of clusters. To be specific, I have refrained to select clusters admixed with the different samples (Figure 6).  </p>
    Finally, I have decided to look into cluster 0, cluster 2, and cluster 12. Cluster 0 seems very specific to COVID-19 patients; Cluster 2 to healthy control; cluster 12 indicates to mostly healthy controls, but moderately admixed with COVID-19 infected subjects. 

**Figure 3: UMAP projection upon groups composed of Healthy Control(HC), Moderate (M), and Severe (S)**
![Figure 3: UMAP projection upon groups composed of Healthy Control(HC), Moderate (M), and Severe (S)](images/analysis/nCoV-umap-group-group.png)

**Figure 4: UMAP projection across samples starting with group names.**
![Figure 4: UMAP projection across samples starting with group names](images/analysis/nCoV-umap-group-sample.png)

**Figure 5: Proportion of clusters upon disease presence.**
![Figure 4: UMAP projection across samples starting with group names](images/analysis/proportion_disease.png)

**Figure 6: UMAP projection with a total of 29 clusters**
![Figure 6: UMAP projection with a total of 29 clusters](images/analysis/umap_markers.png)

### Finding differential markers on the Clusters 0, 2, 12


#### Cluster 0 (COVID19- dominant)
<p> 
    Figure 7 shows dot-plot expression for the most variable genes (Figure 8) across subjects from Severe (0_S) to Healthy Control (0_HC) in cluster 0. To be an instance, WFDC2 gene has been increasingly expressed upon the severeity of disease. 0_NA label most likely correspond to the Moderate (0_M) patients (not shown). There is a list of curated gene by the original authors (i.e. {MARCO, CD48, FCGR3A, TREM2, FCN1, SPP1, FAB4, and CD8A}). Figure 9 shows dot-plot expression for the curated gene lists across subjects. In Cluster 0, FABP4 gene is notable because of the trend of decreasing expression against severeity of disease. 
</p>

#### Cluster 2 (Healthy Control- dominant)
<p> 
    Figure 10 shows dot-plot expression for the most variable genes (Figure 8) across subjects from Severe (0_S) to Healthy Control (0_HC) in cluster 2. A trend of increasing WFDC2 expression against severety of disease has been also observed. Also, dot-plot expression for the curated gene lists has been done (Figure 11). Similar pattern has been observed where FABP4 decreases over the severeity of the disease. </p>


#### Cluster 12 (Admixed) 
<p> 
    Figure 12 shows dot-plot expression for the most variable genes (Figure 8) across subjects from Severe (0_S) to Healthy Control (0_HC) in cluster 12. Those genes are generally expressed exclusively for the COVID-19 patients (0_S, 0_M). No pattern has been observed in the curated gene lists (not shown).  </p> 


**Figure 7: Dot-plot expression in cluster 0**
![Figure 7](images/analysis/cluster0_dotplot.png)


**Figure 8: Most-variable genes**
![Figure 8](images/analysis/cluster0_integrated.png)


**Figure 9: Dot-plot expression in cluster 0 for curated genes**
![Figure 9](images/analysis/cluster0_integrated_curated_dotplot.png)


**Figure 10: Dot-plot expression in cluster 2**
![Figure 10](images/analysis/cluster2_dotplot.png)

**Figure 11: Dot-plot expression in cluster 2 for curated genes**
![Figure 11](images/analysis/cluster2_integrated_curated_dotplot.png)

**Figure 12: Dot-plot expression in cluster 12**
![Figure 12](images/analysis/cluster12_dotplot.png)
