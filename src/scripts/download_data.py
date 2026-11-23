import io
import pandas as pd
import scanpy as sc
import gcsfs

fs = gcsfs.GCSFileSystem()

# GCS bucket path for Virtual Cell Challenge 2025 datasets
vcc_base_path = "gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/"

# load data
train_adata_path = vcc_base_path + "train/adata_Training.h5ad"
val_adata_path = vcc_base_path + "validation/adata_Validation.h5ad"

local_train = "../../data/adata_Training.h5ad"  # make local file
if not os.path.exists(local_train):
    fs.get(train_adata_path, local_train)  # get function downloads the file from google cloud to local path

local_val = "../../data/adata_Validation.h5ad"  # make local file
if not os.path.exists(local_val):
    fs.get(val_adata_path, local_val)  # get function downloads the file from google cloud to local path
