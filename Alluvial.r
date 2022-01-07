library(ggalluvial)
library(dplyr)
library(patchwork)

try<-read.table("try_alluvial.txt",header = T)

mycol<-colorRampPalette(c('#00abef','#64b036','#ffe743','#64b036','#00abef'))(36)

ggplot(try,
       aes(axis1 = tra_v_genes, axis2 = tra_j_genes, axis3 = trb_v_genes, axis4 = trb_j_genes, axis5 = trb_c_genes,
           y= proportion)) +
  scale_x_discrete(limits = c("tra_v_genes", "tra_j_genes", "trb_v_genes", "trb_j_genes", "trb_c_genes"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill=VAC)) +
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=1) +
  theme_void()
  
  
