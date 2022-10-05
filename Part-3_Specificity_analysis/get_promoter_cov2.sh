#!/bin/bash

fantom_dir=/home1/yangyc/share/data_resource/FANTOM

promoter_bed=/home1/GENE_proc/SONGLITING/FANTOM/bed_files/protein_coding.promoter.bed

cd /home1/GENE_proc/SONGLITING/FANTOM/promoter
 
# bedcov
cat ${fantom_dir}/metadata/metadata.FANTOM.dataset.txt |cut  -f5|while read sample

do
	BiosampleID=`echo $sample|cut -d '.' -f3`
	bedtools coverage -s -a $promoter_bed  -b ${fantom_dir}/raw_data/${sample}.hg38.nobarcode.bam > ${BiosampleID}.txt
	#samtools bedcov  $promoter_bed  ${fantom_dir}/raw_data/${sample}.hg38.nobarcode.bam > ${BiosampleID}.txt
	awk -v OFS="\t" '{print $0, "'$BiosampleID'"}' ${BiosampleID}.txt > promoter_${BiosampleID}.txt 
	rm ${BiosampleID}.txt
done
