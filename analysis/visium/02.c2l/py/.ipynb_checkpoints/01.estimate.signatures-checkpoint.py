#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Prepare cell2location reference signatures')
parser.add_argument("infile", type=str,default=None,help='input h5ad file with reference dataset')
parser.add_argument("output", type=str,default=None,help='folder to write output')
parser.add_argument("labels_key", type=str,default=None,help='column in adata.obs to be used as cell type label')
parser.add_argument("--batch_key", type=str,default=None,help='column in adata.obs to be used as bacth (single 10x reaction)')
parser.add_argument("--categorical_covariate_key",default=None, action='append',type=str,help='column in adata.obs to be used as categrical covariates - donor, 3/5, etc (no covariates by default). Multiple columns can be supplied by repetitive usage of this option.')
parser.add_argument("--continuous_covariate_key",default=None, action='append',type=str,help='column in adata.obs to be used as categrical covariates (no covariates by default). Multiple columns can be supplied by repetitive usage of this option.')
parser.add_argument("--gene_id", type=str,default=None,help='column in adata.var to be used as gene id')
parser.add_argument("--cell_count_cutoff", type=int,default=5,help='Gene filtering parameter: All genes detected in less than cell_count_cutoff cells will be excluded.')
parser.add_argument("--cell_percentage_cutoff2", type=float,default=0.03,help='Gene filtering parameter: All genes detected in at least this percentage of cells will be included.')
parser.add_argument("--nonz_mean_cutoff", type=float,default=1.12,help='Gene filtering parameter: genes detected in the number of cells between the above mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff.')
parser.add_argument("--max_epochs", type=int,default=250,help='max_epochs for training')
parser.add_argument("--remove_genes_column", type=str,default=None,help='logical column in adata.var to be used to remove genes, for example mitochonrial. All genes with True in the column will be removed. None (defualt) mean to remove nothing.')

args = parser.parse_args()

import sys
import os
import scanpy as sc
from scipy.sparse import issparse
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import torch
import scvi
from scvi import REGISTRY_KEYS


#######################
# create output folder
os.mkdir(args.output)

sys.stdout = open(args.output+"/c2l.ref.log", "w")
print(args)
print("cuda avaliable: "+str(torch.cuda.is_available()))

# read data
ref = sc.read(args.infile)


mtcnt = np.sum([gene.startswith('MT-') for gene in ref.var.index])
if mtcnt > 0:
    print('There are ' + str(mtcnt) + 'MT genes! Consider to remove them!')

if args.gene_id is not None:
  ref.var[args.gene_id] = ref.var[args.gene_id].astype('string')
  ref.var=ref.var.set_index(args.gene_id)
  print('Raw: cells = '+str(ref.shape[0])+"; genes = " + str(ref.shape[1]))

# filter genes
if args.remove_genes_column != None:
  print('Remove genes by "'+args.remove_genes_column+'". Following genes were removed:')
  print(ref.var[ref.var[args.remove_genes_column]])
  ref = ref[:,~ref.var[args.remove_genes_column]]


# filter genes 
selected = filter_genes(ref,
                        cell_count_cutoff=args.cell_count_cutoff,
                        cell_percentage_cutoff2=args.cell_percentage_cutoff2,
                        nonz_mean_cutoff=args.nonz_mean_cutoff)
                        
plt.savefig(args.output+'/gene.filter.pdf') 

print('Before filtering: cells = '+str(ref.shape[0])+"; genes = " + str(ref.shape[1]))
ref = ref[:, selected].copy()
print('After filtering: cells = '+str(ref.shape[0])+"; genes = " + str(ref.shape[1]))

# train
cell2location.models.RegressionModel.setup_anndata(adata=ref,
                        # 10X reaction / sample / batch
                        batch_key=args.batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=args.labels_key,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=args.categorical_covariate_key,
                        continuous_covariate_keys=args.continuous_covariate_key
                       )

mod = RegressionModel(ref)

mod.view_anndata_setup()

mod.train(max_epochs=args.max_epochs,use_gpu=True,progress_bar_refresh_rate=0)

# plot ELBO loss history during training, removing first 20 epochs from the plot
fig, ax = plt.subplots()
mod.plot_history(20)
plt.savefig(args.output+'/train.history.pdf')

ref = mod.export_posterior(
    ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
mod.save(args.output+"/rsignatures", overwrite=True)
# most likely I do not need this file
ref.write(args.output+"/rsignatures/sc.h5ad")

# save signatures
inf_aver = ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = ref.uns['mod']['factor_names']
inf_aver.to_csv(args.output+'/rsignatures/inf_aver.csv')


# function to plot QCs into file
def plot_QC1(m,plot,summary_name: str = "means",use_n_obs: int = 1000):
  if use_n_obs is not None:
    ind_x = np.random.choice(m.adata_manager.adata.n_obs, np.min((use_n_obs, m.adata.n_obs)), replace=False)
  else:
    ind_x = None
  m.expected_nb_param = m.module.model.compute_expected(
    m.samples[f"post_sample_{summary_name}"], m.adata_manager, ind_x=ind_x
  )
  x_data = m.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]
  if issparse(x_data):
    x_data = np.asarray(x_data.toarray())
  
  mu = m.expected_nb_param["mu"]
  data_node = x_data
  plot.hist2d(np.log10(data_node.flatten()+1), np.log10(mu.flatten()+1), bins=50, norm=mpl.colors.LogNorm())
  plot.set_title("Reconstruction accuracy")
  plot.set(xlabel="Data, log10", ylabel="Posterior sample, values, log10")


def plot_QC2(m,plot,summary_name: str = "means",use_n_obs: int = 1000,scale_average_detection: bool = True):
  inf_aver = m.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
  if scale_average_detection and ("detection_y_c" in list(m.samples[f"post_sample_{summary_name}"].keys())):
    inf_aver = inf_aver * m.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
  aver = m._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
  aver = aver[m.factor_names_]
  plot.hist2d(
    np.log10(aver.values.flatten() + 1),
    np.log10(inf_aver.flatten() + 1),
    bins=50,
    norm=mpl.colors.LogNorm(),)
  plot.set(xlabel="Mean expression for every gene in every cluster", ylabel="Estimated expression for every gene in every cluster")


# unfortunatelly it may not work specifically in case of underpopulated covariates/cell_types. It cannot be fixed on this level, so I'll use "try"
# see https://github.com/BayraktarLab/cell2location/issues/74
fig, (ax1,ax2) = plt.subplots(1,2)
try:
    plot_QC1(mod,plot=ax1,use_n_obs=10000)
except Exception as e:
    print(e)

try:
    plot_QC2(mod,plot=ax2)
except Exception as e:
    print(e)

plt.tight_layout()
plt.savefig(args.output+'/train.QC.pdf')

cell2location.utils.list_imported_modules()
sys.stdout.close()
