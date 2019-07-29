#!/usr/bin/env Rscript

library(ggplot2)
library(SNPRelate)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
tmp=tempfile()
snpgdsVCF2GDS('intermediate_files/vcf_file.QUAL_GT_20.common_snps_only.vcf', 
              tmp,  
              method="biallelic.only")

genofile <- openfn.gds(tmp)

ibs <- snpgdsIBS(genofile, remove.monosnp = TRUE, missing=0.9)
image(ibs$ibs, col=terrain.colors(16))

colnames(ibs$ibs)=ibs$sample.id

ibs_tibble=as.tibble(ibs$ibs)%>%
  mutate(ID1=ibs$sample.id)%>%
  gather('ID2','IBS',c(-ID1))


ibs_aug=ibs_tibble %>%
  mutate(ID1=str_replace(ID1,'-','_'))%>%
  mutate(ID2=str_replace(ID2,'-','_'))


ibs_aug%>%
  ggplot()+
  aes(x=ID1, y=ID2, fill=IBS)+
  geom_raster()+
  scale_fill_viridis_c()+
  theme(axis.text.y = element_text(size=5, family = 'mono'))+
  theme(axis.text.x = element_text(angle=-90, hjust=0, vjust=0, size=5, family='mono'))+
  labs(title='Sequence Similarity Plot',
       x='',
       y='')
print('hello')

ggsave('sequence_similarity.pdf', width = 15,height=15,dpi=300)

write_csv(ibs_aug, 'sequence_similarity.csv')
