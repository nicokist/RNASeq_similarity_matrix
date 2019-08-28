#!/bin/bash
mkdir bam_files
for file in original_bam_files/*.bam; do
	filename=`echo $file | cut -d '/' -f 2 | cut -d "." -f 1`;  
	samtools view -H $file | sed -e 's/SN:chr/SN:/' | samtools reheader - $file > bam_files/${filename}.no_chr.bam; 
done
