library(immunarch)
file_path_TRA = "*"
immdata_mixcr_TRA <- repLoad(file_path)
file_path_TRB = "*"
immdata_mixcr_TRB <- repLoad(file_path)

#Distribution of CDR3 length & Distribution of clonotype abundances & Number of clonotypes
exp_len_TRA <- repExplore(immdata_mixcr_TRA$data, .method = "len", .col = "aa")
exp_cnt_TRA <- repExplore(immdata_mixcr_TRA$data, .method = "count")
exp_vol_TRA <- repExplore(immdata_mixcr_TRA$data, .method = "volume")
exp_len_TRB <- repExplore(immdata_mixcr_TRB$data, .method = "len", .col = "aa")
exp_cnt_TRB <- repExplore(immdata_mixcr_TRB$data, .method = "count")
exp_vol_TRB <- repExplore(immdata_mixcr_TRB$data, .method = "volume")
p_exp_vol_TRA <- vis(exp_vol_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)
p_exp_len_TRA <- vis(exp_len_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)
p_exp_cnt_TRA <- vis(exp_cnt_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)
p_exp_vol_TRB <- vis(exp_vol_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)
p_exp_len_TRB <- vis(exp_len_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)
p_exp_cnt_TRB <- vis(exp_cnt_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)

write.table(exp_len_TRB,file="exp_len_TRB.txt",sep="\t")
(p_exp_vol_TRA+p_exp_vol_TRB)/(p_exp_len_TRA+p_exp_len_TRB)
p_exp_cnt_TRA+p_exp_cnt_TRB

p_exp_vol_TRB+p_exp_len_TRB+p_exp_cnt_TRB


#Clonality
imm_pr_TRA <- repClonality(immdata_mixcr_TRA$data, .method = "clonal.prop")
imm_top_TRA <- repClonality(immdata_mixcr_TRA$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_rare_TRA <- repClonality(immdata_mixcr_TRA$data, .method = "rare")
imm_hom_TRA <- repClonality(immdata_mixcr_TRA$data,
                            .method = "homeo",
                            .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_pr_TRB <- repClonality(immdata_mixcr_TRB$data, .method = "clonal.prop")
imm_top_TRB <- repClonality(immdata_mixcr_TRB$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_rare_TRB <- repClonality(immdata_mixcr_TRB$data, .method = "rare")
imm_hom_TRB <- repClonality(immdata_mixcr_TRB$data,
                            .method = "homeo",
                            .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
vis(imm_top_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)+ vis(imm_top_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)
vis(imm_rare_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)+ vis(imm_rare_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)
vis(imm_hom_TRA, .by = "Status", .meta = immdata_mixcr_TRA$meta)+ vis(imm_hom_TRB, .by = "Status", .meta = immdata_mixcr_TRB$meta)



#Repertoire overlap
imm_ov1 <- repOverlap(immdata_mixcr_TRA$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata_mixcr_TRB$data, .method = "public", .verbose = F)
p8 <- vis(imm_ov1)
p9 <- vis(imm_ov2)
p2 <- vis(imm_ov2, .text.size = 2)

#GU_analysis
geneUsageAnalysis(
  .data=imm_gu_trav_AP, 
  .method = c("js+hclust", "pca+kmeans", "anova", "js+pca+kmeans"),
  .base = 2,
  .norm.entropy = FALSE,
  .cor = c("pearson", "kendall", "spearman"),
  .do.norm = TRUE,
  .laplace = 1e-12,
  .verbose = TRUE,
  .k = 2,
  .eps = 0.01,
  .perp = 1,
  .theta = 0.1
)
#Gene usage
imm_gu_trav<- geneUsage(immdata_mixcr_TRA$data, "hs.trav", .norm = T, .ambig = "exc")
imm_gu_traj <- geneUsage(immdata_mixcr_TRA$data, "hs.traj", .norm = T, .ambig = "exc")
write.table(imm_gu_trav,file="trav_mat.txt",sep="\t")
write.table(imm_gu_traj,file="traj_mat.txt",sep="\t")
  
imm_gu_trbj <- geneUsage(immdata_mixcr_TRB$data, "hs.trbj", .norm = T, .ambig = "exc")
imm_gu_trbv <- geneUsage(immdata_mixcr_TRB$data, "hs.trbv", .norm = T, .ambig = "exc")

write.table(imm_gu_trbv,file="trbv_mat.txt",sep="\t")
write.table(imm_gu_trbj,file="trbj_mat.txt",sep="\t")
##imm_gu_ighv <- geneUsage(immdata_mixcr$data, "hs.ighv", .norm = T)
##fixVis(imm_gu_trav)
##vis(imm_gu_trav, .by = "Status", .meta = immdata_mixcr$meta, .plot = "box")
imm_gu_trav_AP<-imm_gu_trav[c(208,288,627,576,550,391,143,371,541,477,435,71,429,102,587,169,679,569,565,13,340,673,636,270,521,40,315,233,1,460,245,196),]
p_tra_AP<-vis(imm_gu_trav_AP, .by = "Status", .meta = immdata_mixcr_TRA$meta, .plot = "box", .point = FALSE)
p_tra<-vis(imm_gu_trav, .by = "Status", .meta = immdata_mixcr_TRA$meta, .plot = "box", .point = FALSE)
p_trb<-vis(imm_gu_trbv, .by = "Status", .meta = immdata_mixcr_TRB$meta, .plot = "box", .point = FALSE)
vis(imm_gu_trav, .plot = "circos")
fixVis(p_tra)
#vis_box(imm_gu_trav, .by = "Group", .meta = immdata_mixcr$meta, .legend.pos = top)
#Spectratyping
p1 <- vis(spectratype(immdata_mixcr$data[[1]], .quant = "count", .col = "aa+v"))
p2 <- vis(spectratype(immdata_mixcr$data[[2]], .quant = "count", .col = "aa+v"))
p3 <- vis(spectratype(immdata_mixcr$data[[3]], .quant = "count", .col = "aa+v"))
p4 <- vis(spectratype(immdata_mixcr$data[[4]], .quant = "count", .col = "aa+v"))
p1+p2
p3+p4

#kmers
kmers <- getKmers(immdata_mixcr_TRA$data, 5)
vis(kmers)
p1 <- vis(kmers, .head = 5)#top 5
p2 <- vis(kmers, .head = 10)
p3 <- vis(kmers, .head = 30)
(p1 + p2) / p3#Layout with p1 and p2 in a row, p3 in a row
p1 <- vis(kmers, .head = 10, .position = "stack")
p2 <- vis(kmers, .head = 10, .position = "fill")
p3 <- vis(kmers, .head = 10, .position = "dodge")
#Option "stack" stacks all bars on top of each other so you can see the full distribution of kmers. 
#Option "fill" stack all bars on top of each other as well, but normalises it in a such way so you see distribution of counts per-kmer, i.e., you can clearly see which repertoire has more kmer counts than others for a specific kmer. 
#Option "dodge" groups kmer bars of different samples so you can clearly see, which samples has more kmer occurrences overall.
(p1 + p2) / p3#Layout with p1 and p2 in a row, p3 in a row

#Sequence motif visualiztion
kmers_1 <- getKmers(immdata_mixcr$data[[1]], 5)
kp_1 <- kmer_profile(kmers_1, "self")
p5 <- vis(kp)
p1 <- vis(kp_1, .plot = "seq")
kmers_2 <- getKmers(immdata_mixcr$data[[2]], 5)
kp_2 <- kmer_profile(kmers_2, "self")
p6 <- vis(kp)
p2 <- vis(kp_2, .plot = "seq")
kmers_3 <- getKmers(immdata_mixcr$data[[3]], 5)
kp_3<- kmer_profile(kmers_3, "self")
p7 <- vis(kp)
p3 <- vis(kp_3, .plot = "seq")
kmers_4 <- getKmers(immdata_mixcr$data[[4]], 5)
kp_4 <- kmer_profile(kmers_4, "self")
p8 <- vis(kp)
p4 <- vis(kp_4, .plot = "seq")

p1 + p2
