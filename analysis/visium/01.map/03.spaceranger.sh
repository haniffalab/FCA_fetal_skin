#! /bin/bash
#BSUB -G cellgeni
#BSUB -J sr.visium[1-12]
#BSUB -o %J.%I.bsub.log
#BSUB -e %J.%I.bsub.err
#BSUB -q normal
#BSUB -n 16
#BSUB -M62000
#BSUB -R "span[hosts=1] select[mem>62000] rusage[mem=62000]"
#BSUB -q normal

# make command for cellranger-atac
cd /lustre/scratch126/cellgen/cellgeni/pasham/data/2302.fetal.skin/visium/spaceranger

FQSUBDIR=`cat samples.txt | head -n $LSB_JOBINDEX | tail -n 1 | cut -d ' ' -f1`
TAG=`cat samples.txt | head -n $LSB_JOBINDEX | tail -n 1 | cut -d ' ' -f2`

sr=/nfs/cellgeni/pasham/bin/spaceranger-2.0.1/spaceranger
REF=/nfs/cellgeni/pasham/bin/refdata-spaceranger-GRCh38-1.2.0/GRCh38
FQDIR=../fq/$FQSUBDIR
IMG=`ls ../tifs/${TAG}*`  

$sr count \
  --id=$TAG \
  --transcriptome=$REF \
  --fastqs=$FQDIR/$TAG \
  --image=$IMG \
  --sample=$TAG \
  --unknown-slide visium-1\
  --reorient-images true\
  --localcores=16 \
  --localmem=60
