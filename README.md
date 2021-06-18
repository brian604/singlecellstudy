### Singlecellstudy
Project (Soongsil University Graduate Study) and TIL(Today I Learned) on single cell

## COVID project 

- Author: Jiwon Kim 
- For Soongsil University Single Cell sequencing analysis graduate course (under Prof. Junil Kim) 
- Paper/ Data Reference: <br>
[Single-cell landscampe of bronchoalveolar immune cells in patients with COVID-19](https://www.nature.com/articles/s41591-020-0901-9)
- Data that is ready to be analyzed is in data/ in h5 format
- Metadata containing patient (stratification into healthy/moderate/severe) and barcode are in the data/ folder. 
- Folder: image directory is broken down into the followings: 
	- pre_qc : Pre-QCed samples 

#### Introduction of scRNA sequencing analysis
______________

##### Observation of Pre-QC
______________
<p> After downloading the files (h5 format), I have looked at the raw distribution of (1) number of RNA features (2) number of RNA (reads) counts (3) percentage of mitochondria across samples. </p>
I have noticed the followings: <p> 
Normal people (corresponding to the C51, C52, and C100) all have distinct bimodal distribution in RNA feature numbers. On the other hands, COVID-19 infected people have variable distributions. Although it would be nice to do statistical test on the distribution shape (manifold?), I have not done the analysis. <p>
For three moderately infected patients (C141, C142, C144), I have noted that highly diverse distributions of RNA features; from bimodal to long-tailed distribution really close to 0. 
The rest of samples (N = 6) are severely illed patients. Note the general skewed distribution close to 0. Also, RNAs (in terms of both frequencies and kinds) are very diverse, observing from a lot of points concentrated at the bottom for severely illed patients.  </p>

![Figure 1.1: Pre-QC of representative normal subject (C51)](images/pre_qc/C51_qc.png)
![Figure 1.2: Pre-QC of representative moderately illed subject (C142)](images/pre_qc/C142_qc.png)
![Figure 1.3: Pre-QC of representative severely illed subject (C145)](images/pre_qc/C142_qc.png)



 
