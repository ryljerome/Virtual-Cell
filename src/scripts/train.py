# import
import scanpy as sc

# load data
adata = sc.read_h5ad(local_train)  # read the local file into the adata variable
adata_val = sc.read_h5ad(local_val)  # read the local file into the adata variable

# normalize and log transform the data
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize the counts per cell to a total of 10,000 counts per cell
sc.pp.log1p(adata)  # log transform the data which is just ln(data + 1) to avoid taking the log of 0

# create model
# optimizer
# loss

# epoch loop
# train loop

# val loop

### test loop