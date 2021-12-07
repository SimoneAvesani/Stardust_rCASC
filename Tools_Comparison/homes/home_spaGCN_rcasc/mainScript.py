import os,csv,re
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import math
from skimage import io, color
import SpaGCN as spg
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import imagecodecs
import tifffile as tiff 
import anndata
from scanpy import read_10x_h5
from distutils import util


# Read original data and save it to h5ad
print("Read in data")

output_path = sys.argv[1]
use_histo = util.strtobool(sys.argv[2])
use_histo = bool(use_histo)
subsetting = util.strtobool(sys.argv[3])
subsetting = bool(subsetting)
matrix_path = sys.argv[4]
spots_path = sys.argv[5]
image_path = sys.argv[6]
L_resolution = float(sys.argv[7])
n_comps = int(sys.argv[8])
p = sys.argv[9]
if subsetting:
  killed_spots_path = sys.argv[10]

adata = read_10x_h5(matrix_path)
spatial=pd.read_csv(spots_path,sep=",",header=None,na_filter=False,index_col=0)
adata.obs["x1"]=spatial[1]    
adata.obs["x2"]=spatial[2]   
adata.obs["x3"]=spatial[3]   
adata.obs["x4"]=spatial[4]   
adata.obs["x5"]=spatial[5]   

# Select captured samples 
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")

# Read in histology image 
image=io.imread(image_path)
print("Done")

# Read spots for subsampling
if subsetting:
   spots = pd.read_csv(killed_spots_path, header=None) 
else:
  spots = pd.DataFrame(data={"Empty"})

print("Subsampling")
# Transform adata into a dataframe 
adata_df = adata.to_df()

# Subset the dataframe 
adata_df = adata_df[~adata_df.index.isin(spots[0])]
obs = adata.obs[~adata.obs.index.isin(spots[0])]

# Convert dataframe to sparse matrix 
adata_mat = csr_matrix(adata_df)

# Reconstruct annotation data object
adata = anndata.AnnData(X=adata_mat, obs=obs, var=adata.var)

# Set parameters 
b=49
a=1

# Spot coordinates
x2=adata.obs["x2"].tolist()
x3=adata.obs["x3"].tolist()

# Pixel coordinates
x4=adata.obs["x4"].tolist()
x5=adata.obs["x5"].tolist()
adj=spg.calculate_adj_matrix(x=x2,y=x3, x_pixel=x4, y_pixel=x5, beta=b, alpha=a, image=image, histology=use_histo)   


print("Run SpaGCN")
# Set seed
random.seed(200)
torch.manual_seed(200)
np.random.seed(200)

adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=10)
spg.prefilter_specialgenes(adata)

# Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# Set percentage of total expression contributed by neighborhoods
p=float(p)
# Search l 
l=spg.find_l(p=p,adj=adj,start=0.01, end=2,sep=0.005, tol=0.02)

# Define the GCN 
res=0.8
clf=spg.SpaGCN()
clf.set_l(l)

# Train: Init using louvain
clf.train(adata,adj,num_pcs=n_comps,init_spa=True,init="louvain",res=L_resolution,louvain_seed=0,tol=5e-3)

# Predict 
y_pred, prob=clf.predict()

adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')


# Save plot
colors_use=['#111010', '#FFFF00', '#4a6fe3', '#bb7784', '#bec1d4', '#ff9896', '#98df8a', '#ffbb78', '#2ca02c', '#ff7f0e', '#1f77b4', '#800080', '#959595', '#ffff00', '#014d01', '#0000ff', '#ff0000', '#000000']
num_celltype=len(adata.obs["pred"].unique())
adata.uns["pred_colors"]=list(colors_use[:num_celltype])
fig=sc.pl.scatter(adata,alpha=1,x="x5",y="x4",color="pred",show=False,size=120000/adata.shape[0])
fig.set_aspect('equal', 'box')
fig.figure.savefig(output_path + "/Domains.png", dpi=300)
print("Plot saved")

# Save results
clusters = adata.obs["pred"]
clusters = clusters.reset_index()
np.savetxt(output_path + '/Clusters.txt',clusters, delimiter=' ', fmt='%s', header = "Spot Cluster", comments="")
print("Results saved")

exit()
