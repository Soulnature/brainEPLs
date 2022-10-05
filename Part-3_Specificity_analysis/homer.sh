#!bin/bash

# homer: finding motif
# songliting 2020.03.29
# ref:https://mp.weixin.qq.com/s/U7iAS75Mlg6aJFqGp1LfDg

#perl ./motifs/parseJasparMatrix.pl JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt > jaspar.motifs

source ~/miniconda3/bin/activate
#perl ~/miniconda3/bin/findMotifsGenome.pl enhancer_brain.bed  hg38 enhancer_brain  -mknown  /home1/songlt/GTRD/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt

perl ~/miniconda3/bin/findMotifsGenome.pl NonBrain_elevated.bed  hg38 NonBrain_elevated
perl ~/miniconda3/bin/findMotifsGenome.pl Non_specific.bed  hg38 Non_specific
perl ~/miniconda3/bin/findMotifsGenome.pl brain_specific_enhancers.bed  hg38 enhancer_brain_specific
perl ~/miniconda3/bin/findMotifsGenome.pl brain_enhancer.bed  hg38 enhancer_brain 
perl ~/miniconda3/bin/findMotifsGenome.pl enhancer_adult_brain.bed  hg38 enhancer_brain_adult
perl ~/miniconda3/bin/findMotifsGenome.pl enhancer_fetal_brain.bed  hg38 enhancer_brain_fetal
perl ~/miniconda3/bin/findMotifsGenome.pl enhancer_newborn_brain.bed  hg38 enhancer_brain_newborn
