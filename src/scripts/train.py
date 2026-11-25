# import
import numpy as np
import scanpy as sc
import torch
import torch.nn as nn

RUN = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
RANDOM_SEED = 42
TRAIN_VAL_SPLIT = 0.8
HIDDEN_DIM = 512
MAX_EPOCHS = 10
BATCH_SIZE = 256

np.random.seed(RANDOM_SEED)

# load data
adata = sc.read_h5ad('data/adata_Training.h5ad')  # read the local file into adata

# will use the validation data as our testing data
adata_test = sc.read_h5ad('data/adata_Validation.h5ad')

# normalize and log transform the data
sc.pp.normalize_total(adata, target_sum=1e4)  # normalize the counts per cell to a total of 10,000 counts per cell
sc.pp.log1p(adata)  # log transform the data which is just ln(data + 1) to avoid taking the log of 0

# our validation data will be split from the training data
# 80% goes to training and 20% to val
# --> we want to split across perturbations/genes
# otherwise the model would see the same perturbations over and over again
genes = set(adata.obs.target_gene) - {"non-targeting"} # gives all perturbed genes

# sample 80% for training genes first, then take remaining genes for the validation set
train_genes = np.random.choice(list(genes, size=int(len(genes)) * TRAIN_VAL_SPLIT), replace=False)
val_genes = list(genes - set(train_genes))

# also do a sample split to mark which sample is in which group
crtl_idx = adata.obs.target_gene == "non-targeting"
train_idx = adata.obs.target_gene.isin(train_genes)
val_idx = adata.obs.target_gene.isin(val_genes)

X_ctrl = adata[ctrl_idx].X.toarray()  # control expression for the sample
pert_train = adata[train_idx].obs.target_gene.values  # list of the training samples (the 80%)
y_train = adata[train_idx].X.toarray()

pert_val = adata[val_idx].obs.target_gene.values
y_val = adata[val_idx].X.toarray()

## need to convert our data into numbers
# our expression data is numerical but the gene names are strings
# so how do encode the names as numbers so the model can be trained on one set 
# of genes, and then validated on a new set of genes
# --> use embeddings: represent each gene as a vector of a fixed length
# the vector reflects relationships between the genes

# create model
model_config = {
    "expr_dim":X_ctrl.shape[1],
    "emb_dim":...,
    "hidden_dim":HIDDEN_DIM,
}

expr_dim = X_ctrl.shape[1]
emb_dim = ...

# architecture
# expr + emb_dim -> hidden_dim -> expr_dim
#model = VirtualCellModel(expr_dim=expr_dim, emb_dim=emb_dim, hidden_dim=HIDDEN_DIM)
model = VirtualCellModel(**model_config)

# optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

# loss
criterion = nn.L1Loss()

# epoch loop
for epoch in range(MAX_EPOCHS):
    # train loop
    train_loss = 0.0
    model.train()
    for batch_ctrl_expr, batch_pert, batch_target_expr in train_loader:
        optimizer.zero_grad()

        # providing the model the batch control expression and perturbation, 
        # should give us back predicted expression
        pred_expr = model(batch_ctrl_expr, batch_pert)

        loss = criterion(pred_expr, batch_target_expr)
        train_loss += loss.item() * batch_ctrl_expr.size(0) # running sum

        loss.backward()
        optimizer.step()

    train_loss /= len(train_loader.dataset)

    # val loop
    val_loss = 0.0
    model.eval()
    with torch.no_grad():
        for batch_ctrl_expr, batch_pert, batch_target_expr in val_loader:

            pred_expr = model(batch_ctrl_expr, batch_pert)
            
            loss = criterion(pred_expr, batch_target_expr)
            val_loss += loss.item() * batch_ctrl_expr.size(0)

    val_loss /= len(val_loader.dataset)

    print(f"Epoch {epoch+1}/{MAX_EPOCHS} - Train Loss: {train_loss:.4f} - Val Loss: {val_loss:.4f}")

    # save model checkpoint
    torch.save({
        "model_state_dict":model.state_dict(),
        "model_config":model_config}, 
        f"model_epoch_{RUN}.pth")

### test loop