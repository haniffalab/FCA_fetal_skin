# Aim

Describe the various analyses performed.

**To be completed**

# Trajectory analysis

Trajectory analysis was performed with Monocle3 for three broad cell types separately:

* keratinocytes, see kc_trajectory.R or variant **TODO: check**  
* mesenchymal, see ms_trajectory_20220322.R
* vascular endothelium, see ve_trajectory_20220503.R

## vascular endothelium

**To be completed**

See ve_trajectory_20220503.R,
which Ni first wrote to analyse the fetal skin data set 'fsk' then the 'pooled' data set.

The input file ve_trajectory_20220503.R read in is
pooled.vascular_endothelium.count_with_PCA_UMAP_for_monocle.20220531.h5ad,
which would be made by prepare_endo_monocle_input.ipynb,
which itself reads in pooled_endothelium.processed.h5ad,
which would be created by 2_melanocytes.ipynb (May 20  2020).
There is a 5_endothelium.ipynb (May 27  2022) but that checks pooled_endothelium.processed.h5ad
but does not create it.

### Create pooled_endothelium.processed.h5ad 

* code: 2_melanocytes.ipynb
* aim: merge fetal and organoid data sets for separate broad cell types:
keratinocytes, melanocytes, endothelium, neuronal, 
* input:
  * organoid.cellxgene.h5ad
  * 20200403_post_annot3_cleanup/fetal_skin.keratinocytes.doublet_removed_processed.20200403.h5ad
  * 20200403_post_annot3_cleanup/fetal_skin.melanocytes.doublet_removed_processed.20200403.h5ad
  * 20200403_post_annot3_cleanup/fetal_skin.endothelium.doublet_removed_processed.20200403.h5ad
  * data/h5ad/20200114/fetal_skin.stroma.doublet_removed_processed.20200114.h5ad (for neuronal and stroma)
* note: uses `scanpy_scripts` now replaced by `sctk`.
* output:
  * pooled_keratinocytes.processed.h5ad
  * pooled_melanocytes.processed.h5ad
  * pooled_endothelium.processed.h5ad
  * pooled_neuronal.processed.h5ad
  * pooled_stroma.processed.h5ad

<!--
# Trajectory analysis
-->

