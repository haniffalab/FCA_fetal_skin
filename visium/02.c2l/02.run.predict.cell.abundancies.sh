#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J c2l.pred[1]
#BSUB -o %J.%I.c2l.pred.log
#BSUB -e %J.%I.c2l.pred.err
#BSUB -n 2
#BSUB -M64000
#BSUB -R "span[hosts=1] select[mem>64000] rusage[mem=64000]"

# for gpu-normal
#BSUB -q gpu-normal
#BSUB -gpu "mode=shared:j_exclusive=yes:gmem=22000:num=1"

# for gpu-cellgeni-a100
##BSUB -q gpu-cellgeni-a100
##BSUB -m dgx-b11
##BSUB -gpu "mode=shared:j_exclusive=no:gmem=62000:num=1"

# uncomment (and edit paths) if you run it not as cellgeni-su
export SINGULARITY_CACHEDIR=/nfs/cellgeni/pasham/singularity/cache
export SINGULARITY_PULLFOLDER=/nfs/cellgeni/pasham/singularity
export SINGULARITY_TMPDIR=/nfs/cellgeni/pasham/singularity/tmp
export PATH=/software/singularity-v3.6.4/bin:$PATH


cd /nfs/cellgeni/pasham/projects/2302.fetal.skin/data.nfs/my/visium/skin.c2l
WDIR=`pwd -P`
IMAGE=/nfs/cellgeni/singularity/images/c2l.jhub.221206.v0.1.sif # based on commit 36e4f007e8fba4cb85c13b9bff47a4f6fbae9295

c2lpred=/nfs/cellgeni/pasham/projects/2302.fetal.skin/src/my/visium/02.c2l/py/02.predict.cell.abundancies.py

# edit here
# by default script will submit array job (edit -J bsub option above) for all combinations of refs and alphas specified below
refs=(fetal_skin.norm.maternal_removed.20220202.filtered10.denorm fetal_skin.norm.maternal_removed.20220202.filtered.denorm)
alphas=(20)
i=$(($LSB_JOBINDEX-1))

alp=${alphas[$((i / ${#refs[@]}))]}
ref=${refs[$((i % ${#refs[@]}))]}

# script expects reference signaatures to be stored in ref subfolder, but paths can be specified directly here
REF=ref/${ref}/rsignatures/inf_aver.csv
OUT=pred/${ref}.${alp}
VISIN=viss.h5ad # path to combined visum h5ad


# edit below if you want to change some defaults (ncells, epochs)
singularity exec --nv --bind /lustre,/nfs $IMAGE /bin/bash -c "nvidia-smi;cd ${WDIR}; \
 $c2lpred \
  --batch_key library_id \
  --detection_alpha ${alp} \
  --N_cells_per_location 30 \
  --max_epochs 50000 \
  $VISIN \
  $REF \
  $OUT"
