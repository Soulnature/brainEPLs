#!/bin/bash

fantom_dir=/home1/yangyc/share/data_resource/FANTOM

enhancer_bed=/home1/GENE_proc/SONGLITING/FANTOM/bed_files/ELS.bed

cd /home1/GENE_proc/SONGLITING/FANTOM/enhancer
 
# bedcov
cat ${fantom_dir}/metadata/metadata.FANTOM.dataset.txt |sed -n '2,$p'|cut  -f5|while read sample

#cat todo.txt |cut  -f5|while read sample
do
	BiosampleID=`echo $sample|cut -d '.' -f3`
	bedtools coverage -a $enhancer_bed  -b ${fantom_dir}/raw_data/${sample}.hg38.nobarcode.bam > ${BiosampleID}.txt
	awk -v OFS="\t" '{print $0, "'$BiosampleID'"}' ${BiosampleID}.txt > enhancer_${BiosampleID}.txt 
	rm ${BiosampleID}.txt
done
