import stlearn as st
from pathlib import Path
import sys
from distutils import util
import random
import string
import numpy as np
import pandas as pd

print("Read in data")
input_path = sys.argv[1]
output_path = sys.argv[2]
n_comps = int(sys.argv[3])
L_resolution = float(sys.argv[4])
subsetting = util.strtobool(sys.argv[5])
subsetting = bool(subsetting)
if subsetting:
  killed_spots_path = sys.argv[6]

letters = string.ascii_lowercase
tile_random_path = ''.join(random.choice(letters) for i in range(5))

TILE_PATH = Path("/tmp/" + tile_random_path)
TILE_PATH.mkdir(parents=True, exist_ok=True)

data = st.Read10X(input_path)
st.pp.filter_genes(data,min_cells=10)
st.pp.normalize_total(data)
st.pp.log1p(data)
st.pp.tiling(data, TILE_PATH)
st.pp.extract_feature(data)
st.em.run_pca(data,n_comps=n_comps)
st.spatial.SME.SME_normalize(data, use_data="raw")
data.X = data.obsm['raw_SME_normalized']

if subsetting:
   spots = pd.read_csv(killed_spots_path, header=None) 
else:
  spots = pd.DataFrame(data={"Empty"})

keep = set(data.obs.index.values).difference(set(spots[0]))
keep = list(keep)

data = data[keep,]
st.pp.scale(data)
st.em.run_pca(data,n_comps=n_comps)
# n_neighbors = 15 is the default value
st.pp.neighbors(data,n_neighbors=15,use_rep='X_pca')
st.tl.clustering.louvain(data, resolution=L_resolution)

clusters = data.obs['louvain']
clusters = clusters.sort_index()
clusters.to_csv(output_path + '/Clusters.txt',sep="\t")
