#! /bin/bash -e
#BSUB -G cellgeni
#BSUB -J c2l.ref[1]
#BSUB -o %J.%I.c2l.ref.log
#BSUB -e %J.%I.c2l.ref.err
#BSUB -n 2
#BSUB -M64000
#BSUB -R "span[hosts=1] select[mem>64000] rusage[mem=64000]"

#BSUB -q gpu-normal
#BSUB -gpu "mode=shared:j_exclusive=no:gmem=22000:num=1"

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
c2lref=/nfs/cellgeni/pasham/projects/2302.fetal.skin/src/my/visium/02.c2l/py/01.estimate.signatures.py

# edit below
# script expects all ref h5ads to be in working directory. It will submit array job (edit -J bsub option above), one instance per ref h5ad
# list here names of ref  files (without ".h5ad"). Estimated signatures will be named in the same way
ref=(fetal_skin.norm.maternal_removed.20220202.filtered10.denorm fetal_skin.norm.maternal_removed.20220202.filtered.denorm)

# set parameters
REFOUT=ref/${ref[$LSB_JOBINDEX-1]}
REFIN=${ref[$LSB_JOBINDEX-1]}.h5ad

CELLTYPE="joint_annotation_20220202" # names of adata.obs column with celltype annotation
SAMPLE="sanger_id" # names of adata.obs column with information about 10x sample (to be used as batch in c2l)

# usual set of parameters is defined above, but you may need to edit code below as well. See manual in py/01.estimate.signatures.py
singularity exec --nv --bind /lustre,/nfs $IMAGE /bin/bash -c "nvidia-smi;cd ${WDIR}; \
 $c2lref \
  --batch_key $SAMPLE \
  --categorical_covariate_key donor \
  --categorical_covariate_key chemistry \
  $REFIN \
  $REFOUT \
  $CELLTYPE"
