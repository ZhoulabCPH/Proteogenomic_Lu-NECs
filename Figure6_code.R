###   Disease progression-associated proteogenomic alterations for oncogenic pathways  ###

###  Figure 6A  ###
library(pheatmap)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
annotation_col <- data.frame(type = clini$type,
                             ajcc = as.factor(clini$AJCCstage),
                             T_stage = as.factor(clini$T_stage),
                             N_metastasis = as.factor(clini$N_metastasis),
                             M_metastasis = as.factor(clini$M_metastasis),
                             Postoperative_adjuvant_chemotherapy = as.factor(clini$Postoperative_adjuvant_chemotherapy))
rownames(annotation_col) <- rownames(clini)
ann_colors = list(type = c("LCNEC" = "#526187","LCNEC+SCLC" = "#ED8A3F",
                           "SCLC" = "#9B8281"),
                  ajcc = c("I" = "#F7F2E9","II" = "#D2D2D2",
                           "III" = "#A9A9A9","IV" = "#848484"),
                  T_stage = c("1"="#dee2d1","2"="#91c0c0",
                              "3"="#c29e2e","4"="#647370"),
                  N_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  M_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy = c("No"="#A1DFDB","Yes"="#F79990"))
linshi <- t(data.frame(x = 1:93, y = 1:93))
colnames(linshi) <- rownames(clini)
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
pheatmap(linshi,fontsize=6,cellheight = 1,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = T)



###   Figure 6B   ###
library(survival)
library(survminer)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical <- clinical[which(clinical$Postoperative_adjuvant_chemotherapy == "Yes"),]
# AJCC stage
all.aj.dfs <- survfit(Surv(DFSm,DFSstate)~AJCCstage,
                     data = clinical)
summary(all.aj.dfs)
ggsurvplot(all.aj.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#F7F2E9","#D2D2D2",
                       "#A9A9A9","#848484"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("I","II","III","IV"), 
           xlab = "Time (months)",ylab = "Disease-free survival")
# N_metastasis
all.n.dfs <- survfit(Surv(DFSm,DFSstate)~N_metastasis,
                     data = clinical)
summary(all.n.dfs)
ggsurvplot(all.n.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#D3E2F2","#F79990"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("No","Yes"), 
           xlab = "Time (months)",ylab = "Disease-free survival")


###  Figure 6C   ###
library(ggplot2)
library(ggpubr)
library(reshape2)
library(gghalves)
library(clinfun)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
inter_name <- intersect(colnames(data.origin),rownames(clinical))
clinical <- clinical[inter_name,]
ajcc1 <- clinical[which(clinical$AJCCstage == "I"),]
ajcc2 <- clinical[which(clinical$AJCCstage == "II"),]
ajcc3 <- clinical[which(clinical$AJCCstage == "III"),]
ajcc4 <- clinical[which(clinical$AJCCstage == "IV"),]
ajcc1.PL <- ajcc1[which(substr(rownames(ajcc1),1,2) == "PL"),]
ajcc2.PL <- ajcc2[which(substr(rownames(ajcc2),1,2) == "PL"),]
ajcc3.PL <- ajcc3[which(substr(rownames(ajcc3),1,2) == "PL"),]
ajcc1.LS <- ajcc1[which(substr(rownames(ajcc1),1,2) == "LS"),]
ajcc2.LS <- ajcc2[which(substr(rownames(ajcc2),1,2) == "LS"),]
ajcc3.LS <- ajcc3[which(substr(rownames(ajcc3),1,2) == "LS"),]
ajcc1.PS <- ajcc1[which(substr(rownames(ajcc1),1,2) == "PS"),]
ajcc2.PS <- ajcc2[which(substr(rownames(ajcc2),1,2) == "PS"),]
ajcc3.PS <- ajcc3[which(substr(rownames(ajcc3),1,2) == "PS"),]
data.PL <- data.frame(i = apply(data.origin[,rownames(ajcc1.PL)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.PL)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.PL)],1,median))
data.PL <- as.data.frame(t(data.PL))
data.LS <- data.frame(i = apply(data.origin[,rownames(ajcc1.LS)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.LS)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.LS)],1,median))
data.LS <- as.data.frame(t(data.LS))
data.PS <- data.frame(i = apply(data.origin[,rownames(ajcc1.PS)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.PS)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.PS)],1,median))
data.PS <- as.data.frame(t(data.PS))
result.PL <- read.csv("AJCC_LCNEC.csv",stringsAsFactors = F,
                   row.names = 1,check.names = F)
result.LS <- read.csv("AJCC_LS.csv",stringsAsFactors = F,
                      row.names = 1,check.names = F)
result.PS <- read.csv("AJCC_SCLC.csv",stringsAsFactors = F,
                      row.names = 1,check.names = F)
## ALL merge
data.PL <- data.frame(i = apply(data.origin[,rownames(ajcc1.PL)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.PL)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.PL)],1,median))
data.PL <- data.PL[rownames(result.PL),]
data.PL[rownames(result.PL[which(result.PL$trend != "no"),]),]
data.PL$trend <- result.PL$trend
data.PL$protein <- rownames(data.PL)
data.PL$cancer_type <- rep("PL",dim(data.PL)[1])
data.LS <- data.frame(i = apply(data.origin[,rownames(ajcc1.LS)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.LS)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.LS)],1,median))
data.LS <- data.LS[rownames(result.LS),]
data.LS[rownames(result.LS[which(result.LS$trend != "no"),]),]
data.LS$trend <- result.LS$trend
data.LS$protein <- rownames(data.LS)
data.LS$cancer_type <- rep("LS",dim(data.LS)[1])
data.PS <- data.frame(i = apply(data.origin[,rownames(ajcc1.PS)],1,median),
                      ii = apply(data.origin[,rownames(ajcc2.PS)],1,median),
                      iii = apply(data.origin[,rownames(ajcc3.PS)],1,median))
data.PS <- data.PS[rownames(result.PS),]
data.PS[rownames(result.PS[which(result.PS$trend != "no"),]),]
data.PS$trend <- result.PS$trend
data.PS$protein <- rownames(data.PS)
data.PS$cancer_type <- rep("PS",dim(data.PS)[1])

data.merge <- rbind(data.PL[which(data.PL$trend != "no"),],
                    data.LS[which(data.LS$trend != "no"),],
                    data.PS[which(data.PS$trend != "no"),])
dim(data.merge)
length(unique(rownames(data.merge)))
data.merge1 <- reshape2::melt(data.merge,id = c('protein',"trend","cancer_type"))
data.merge1$cancer_type <- factor(data.merge1$cancer_type,
                                  levels = c("PL","LS","PS"))
data.merge1.up <- data.merge1[which(data.merge1$trend == "dizeng"),]
data.merge1.do <- data.merge1[which(data.merge1$trend == "dijian"),]
data.merge1.up.PL <- data.merge1.up[which(data.merge1.up$cancer_type == "PL"),]
data.merge1.up.LS <- data.merge1.up[which(data.merge1.up$cancer_type == "LS"),]
data.merge1.up.PS <- data.merge1.up[which(data.merge1.up$cancer_type == "PS"),]
data.merge1.do.PL <- data.merge1.do[which(data.merge1.do$cancer_type == "PL"),]
data.merge1.do.LS <- data.merge1.do[which(data.merge1.do$cancer_type == "LS"),]
data.merge1.do.PS <- data.merge1.do[which(data.merge1.do$cancer_type == "PS"),]
## LCNEC
ggplot(data = data.merge1.up.PL,aes(x = variable,y = value,
                          group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#8A9EB5')) +
  theme_bw()
ggplot(data = data.merge1.do.PL,aes(x = variable,y = value,
                                 group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#8A9EB5')) +
  theme_bw()
## C-SCLC
ggplot(data = data.merge1.up.LS,aes(x = variable,y = value,
                                    group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#ED8A3F')) +
  theme_bw()
ggplot(data = data.merge1.do.LS,aes(x = variable,y = value,
                                    group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#ED8A3F')) +
  theme_bw()
## SCLC
ggplot(data = data.merge1.up.PS,aes(x = variable,y = value,
                                    group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#9B8281')) +
  theme_bw()
ggplot(data = data.merge1.do.PS,aes(x = variable,y = value,
                                    group = protein,color = cancer_type)) +
  geom_line()+ 
  scale_color_manual(values = c('#9B8281')) +
  theme_bw()



###  Figure 6C  heatmap  ###
library(pheatmap)
up <- read.csv("up.csv",
                   stringsAsFactors = F,row.names = 1,check.names = F)
down <- read.csv("down.csv",
                   stringsAsFactors = F,row.names = 1,check.names = F)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical <- clinical[which(clinical$AJCCstage != "IV"),]
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
inter_name <- intersect(colnames(data.origin),rownames(clinical))
data.origin <- data.origin[,inter_name]
clinical <- clinical[inter_name,]
clinical.pl <- clinical[which(substr(rownames(clinical),1,2) == "PL"),]
clinical.ls <- clinical[which(substr(rownames(clinical),1,2) == "LS"),]
clinical.ps <- clinical[which(substr(rownames(clinical),1,2) == "PS"),]

clinical.pl <- clinical.pl[order(clinical.pl$AJCCstage),]
clinical.ls <- clinical.ls[order(clinical.ls$AJCCstage),]
clinical.ps <- clinical.ps[order(clinical.ps$AJCCstage),]

data.origin.pl <- data.origin[c(rownames(dijian[which(dijian$cancer_type == "PL"),]),
                                rownames(dizeng[which(dizeng$cancer_type == "PL"),])),
                              rownames(clinical.pl)]
data.origin.ls <- data.origin[c(rownames(dijian[which(dijian$cancer_type == "LS"),]),
                                rownames(dizeng[which(dizeng$cancer_type == "LS"),])),
                              rownames(clinical.ls)]
data.origin.ps <- data.origin[c(rownames(dijian[which(dijian$cancer_type == "PS"),]),
                                rownames(dizeng[which(dizeng$cancer_type == "PS"),])),
                              rownames(clinical.ps)]
linshi.pl <- apply(data.origin.pl, 1, scale)
linshi.pl <- t(linshi.pl)
colnames(linshi.pl) <- colnames(data.origin.pl)
linshi.pl[linshi.pl < (-2)] <- (-2)
linshi.pl[linshi.pl > 2] <- 2
linshi.ls <- apply(data.origin.ls, 1, scale)
linshi.ls <- t(linshi.ls)
colnames(linshi.ls) <- colnames(data.origin.ls)
linshi.ls[linshi.ls < (-2)] <- (-2)
linshi.ls[linshi.ls > 2] <- 2
linshi.ps <- apply(data.origin.ps, 1, scale)
linshi.ps <- t(linshi.ps)
colnames(linshi.ps) <- colnames(data.origin.ps)
linshi.ps[linshi.ps < (-2)] <- (-2)
linshi.ps[linshi.ps > 2] <- 2
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
annotation_col <- data.frame(ajcc = clinical$AJCCstage,
                             N_metastasis = as.factor(clinical$N_metastasis))
rownames(annotation_col) <- rownames(clinical)
ann_colors = list(N_metastasis = c("No"="#D3E2F2","Yes"="#F79990"),
                  ajcc = c("I" = "#F7F2E9","II" = "#D2D2D2",
                           "III" = "#A9A9A9"))
zengjian <- rbind(dizeng,dijian)
rownames(linshi.pl) <- zengjian[rownames(linshi.pl),7]
rownames(linshi.ls) <- zengjian[rownames(linshi.ls),7]
rownames(linshi.ps) <- zengjian[rownames(linshi.ps),7]
pheatmap(linshi.pl,fontsize=6,cellwidth = 3,cellheight = 9,
         color  = colorRampPalette(c(blue,white,red))(100),
         clustering_method = "ward.D2",
         border_color = "grey60",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
pheatmap(linshi.ls,fontsize=6,cellwidth = 3,cellheight = 9,
         color  = colorRampPalette(c(blue,white,red))(100),
         clustering_method = "ward.D2",
         border_color = "grey60",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
pheatmap(linshi.ps,fontsize=6,cellwidth = 3,cellheight = 9,
         color  = colorRampPalette(c(blue,white,red))(100),
         clustering_method = "ward.D2",
         border_color = "grey60",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
 

###   Figure 6D     ###
library(ggplot2)
library(ggpubr)
library(gghalves)
library(ggbreak)
up <- read.csv("up.csv",
                   stringsAsFactors = F,row.names = 1,check.names = F)
down <- read.csv("down.csv",
                   stringsAsFactors = F,row.names = 1,check.names = F)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical <- clinical[which(clinical$AJCCstage != "IV"),]
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
inter_name <- intersect(colnames(data.origin),rownames(clinical))
data.origin <- data.origin[,inter_name]
clinical <- clinical[inter_name,]
clinical.pl <- clinical[which(substr(rownames(clinical),1,2) == "PL"),]
clinical.ls <- clinical[which(substr(rownames(clinical),1,2) == "LS"),]
clinical.ps <- clinical[which(substr(rownames(clinical),1,2) == "PS"),]

clinical.pl <- clinical.pl[order(clinical.pl$AJCCstage),]
clinical.ls <- clinical.ls[order(clinical.ls$AJCCstage),]
clinical.ps <- clinical.ps[order(clinical.ps$AJCCstage),]

data.origin.pl <- data.origin[c(rownames(dizeng[which(dizeng$cancer_type == "PL"),])),
                              rownames(clinical.pl)]
data.origin.ls <- data.origin[c(rownames(dizeng[which(dizeng$cancer_type == "LS"),])),
                              rownames(clinical.ls)]
data.origin.ps <- data.origin[c(rownames(dizeng[which(dizeng$cancer_type == "PS"),])),
                              rownames(clinical.ps)]
data.origin.pl <- as.data.frame(t(data.origin.pl))
data.origin.ls <- as.data.frame(t(data.origin.ls))
data.origin.ps <- as.data.frame(t(data.origin.ps))
data.origin.pl$n_stage <- clinical.pl$N_metastasis
data.origin.ls$n_stage <- clinical.ls$N_metastasis
data.origin.ps$n_stage <- clinical.ps$N_metastasis
wilcox.test(data.origin.pl[,"Q9BXX0"]~data.origin.pl$n_stage)
wilcox.test(data.origin.ls[,"O95989"]~data.origin.ls$n_stage)
wilcox.test(data.origin.ls[,"O00746"]~data.origin.ls$n_stage)
wilcox.test(data.origin.ps[,"P10645"]~data.origin.ps$n_stage)
wilcox.test(data.origin.ps[,"P19022"]~data.origin.ps$n_stage)
wilcox.test(data.origin.ps[,"P61764"]~data.origin.ps$n_stage)

data.box <- data.frame(y = c(data.origin.pl$Q9BXX0,
                             data.origin.ls$O00746,
                             data.origin.ls$O95989,
                             data.origin.ps$P10645,
                             data.origin.ps$P19022,
                             data.origin.ps$P61764),
                       x = c(rep("Q9BXX0",dim(data.origin.pl)[1]),
                             rep("O00746",dim(data.origin.ls)[1]),
                             rep("O95989",dim(data.origin.ls)[1]),
                             rep("P10645",dim(data.origin.ps)[1]),
                             rep("P19022",dim(data.origin.ps)[1]),
                             rep("P61764",dim(data.origin.ps)[1])),
                       n_stage = c(data.origin.pl$n_stage,
                                   data.origin.ls$n_stage,
                                   data.origin.ls$n_stage,
                                   data.origin.ps$n_stage,
                                   data.origin.ps$n_stage,
                                   data.origin.ps$n_stage))
ggplot() +
  geom_half_boxplot(data = data.box, 
                    aes(x=x, y=y, 
                        fill = n_stage),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data.box, 
                  aes(x=x, y=y, color=n_stage), 
                  shape = 16, alpha=1, size = 0.8,
                  transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#D3E2F2','#F79990'))+
  scale_color_manual(values = c('#D3E2F2','#F79990'))+
  stat_compare_means(data = data.box, 
                     aes(x=x, y=y),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")




###  Figure 6E  ###
library(ggplot2)
library(Hmisc)
library(ambient)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical <- clinical[which(clinical$AJCCstage != "IV"),]
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
data.pl <- data.origin["Q9BXX0",
                       which(substr(colnames(data.origin),1,2) == "PL")]
data.ls <- data.origin[c("O95989","O00746"),
                       which(substr(colnames(data.origin),1,2) == "LS")]
data.ps <- data.origin[c("P61764","P19022","P10645"),
                       which(substr(colnames(data.origin),1,2) == "PS")]
data.pl <- as.data.frame(t(data.pl))
data.ls <- as.data.frame(t(data.ls))
data.ps <- as.data.frame(t(data.ps))
fge.score <- read.table("29_FGE.txt",sep = "\t",check.names = F,
                        stringsAsFactors = F,row.names = 1)
fge.score <- as.data.frame(t(fge.score))
fge.score.pl <- fge.score[rownames(data.pl),]
fge.score.ls <- fge.score[rownames(data.ls),]
fge.score.ps <- fge.score[rownames(data.ps),]
### PL
data.pl$emt <- fge.score.pl$`EMT signature`
data.pl$tpr <- fge.score.pl$`Tumor proliferation rate`
data.pl$caf <- fge.score.pl$`Cancer-associated fibroblasts`
data.pl$matrix <- fge.score.pl$Matrix
data.pl$matrix.remo <- fge.score.pl$`Matrix remodeling`
pl.cor <- rcorr(as.matrix(data.pl),type = "spearman")
pl.cor.p <- pl.cor$P
pl.cor.r <- pl.cor$r
clinical.pl <- clinical[rownames(data.pl),]
data.pl$ajcc <- clinical.pl$AJCCstage
### LS
data.ls$emt <- fge.score.ls$`EMT signature`
data.ls$tpr <- fge.score.ls$`Tumor proliferation rate`
data.ls$caf <- fge.score.ls$`Cancer-associated fibroblasts`
data.ls$matrix <- fge.score.ls$Matrix
data.ls$matrix.remo <- fge.score.ls$`Matrix remodeling`
ls.cor <- rcorr(as.matrix(data.ls),type = "spearman")
ls.cor.p <- ls.cor$P
ls.cor.r <- ls.cor$r
clinical.ls <- clinical[rownames(data.ls),]
data.ls$ajcc <- clinical.ls$AJCCstage
### PS
data.ps$emt <- fge.score.ps$`EMT signature`
data.ps$tpr <- fge.score.ps$`Tumor proliferation rate`
data.ps$caf <- fge.score.ps$`Cancer-associated fibroblasts`
data.ps$matrix <- fge.score.ps$Matrix
data.ps$matrix.remo <- fge.score.ps$`Matrix remodeling`
ps.cor <- rcorr(as.matrix(data.ps),type = "spearman")
ps.cor.p <- ps.cor$P
ps.cor.r <- ps.cor$r
clinical.ps <- clinical[rownames(data.ps),]
data.ps$ajcc <- clinical.ps$AJCCstage
pl.cor.r <- as.data.frame(pl.cor.r)
ls.cor.r <- as.data.frame(ls.cor.r)
ps.cor.r <- as.data.frame(ps.cor.r)
pl.cor.p <- as.data.frame(pl.cor.p)
ls.cor.p <- as.data.frame(ls.cor.p)
ps.cor.p <- as.data.frame(ps.cor.p)
library(circlize)
library(ggplot2)
library(ggforce)
library(ComplexHeatmap)
cor.r.merge <- cbind(pl.cor.r[2:6,1],ls.cor.r[3:7,1:2],ps.cor.r[4:8,1:3])
cor.p.merge <- cbind(pl.cor.p[2:6,1],ls.cor.p[3:7,1:2],ps.cor.p[4:8,1:3])
colnames(cor.r.merge)[1] <- "Q9BXX0"
colnames(cor.p.merge)[1] <- "Q9BXX0"
cor.r.merge <- t(cor.r.merge)
cor.p.merge <- t(cor.p.merge)
col_fun <- colorRamp2(c(-1,0, 1), 
                      c("#1F2A38","#F3F3F3","#FFCD02"))
circos.clear()
circos.par(gap.after = c(90))
circos.heatmap(cor.r.merge, 
               rownames.side = "outside",
               col = col_fun,track.height = 0.7)
grid.draw(Legend(title = "Title", col_fun = col_fun))
circos.clear()


###  Figure 6F   ###
library(maxstat)
library(survival)
library(survminer)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical <- clinical[which(clinical$AJCCstage != "IV"),]
clinical <- clinical[which(clinical$Postoperative_adjuvant_chemotherapy == "Yes"),]
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
data.pl <- data.origin["Q9BXX0",
                       which(substr(colnames(data.origin),1,2) == "PL")]
data.ps <- data.origin[c("P61764","P19022"),
                       which(substr(colnames(data.origin),1,2) == "PS")]
clinical.pl <- clinical[colnames(data.pl),]
clinical.ps <- clinical[colnames(data.ps),]
data.pl <- as.data.frame(t(data.pl))
data.ps <- as.data.frame(t(data.ps))
clinical.pl <- cbind(clinical.pl,data.pl)
clinical.ps <- cbind(clinical.ps,data.ps)
EMILIN2 <- maxstat.test(Surv(DFSm,DFSstate) ~ Q9BXX0,
                        data = clinical.pl, smethod = "LogRank")
STXBP1 <- maxstat.test(Surv(DFSm,DFSstate) ~ P61764,
                       data = clinical.ps, smethod = "LogRank")
CDH2 <- maxstat.test(Surv(DFSm,DFSstate) ~ P19022, 
                     data = clinical.ps, smethod = "LogRank")
clinical.pl$EMILIN2.lab <- clinical.pl$Q9BXX0
clinical.pl$EMILIN2.lab[clinical.pl$Q9BXX0 <= EMILIN2$estimate] <- "Low"
clinical.pl$EMILIN2.lab[clinical.pl$Q9BXX0 > EMILIN2$estimate] <- "High"
clinical.ps$STXBP1.lab <- clinical.ps$P61764
clinical.ps$STXBP1.lab[clinical.ps$P61764 <= STXBP1$estimate] <- "Low"
clinical.ps$STXBP1.lab[clinical.ps$P61764 > STXBP1$estimate] <- "High"
clinical.ps$CDH2.lab <- clinical.ps$P19022
clinical.ps$CDH2.lab[clinical.ps$P19022 <= CDH2$estimate] <- "Low"
clinical.ps$CDH2.lab[clinical.ps$P19022 > CDH2$estimate] <- "High"
EMILIN2.dfs <- survfit(Surv(DFSm,DFSstate)~EMILIN2.lab,
                   data = clinical.pl)
summary(EMILIN2.dfs)
ggsurvplot(EMILIN2.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#FF6361","#CFD3CE"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("High","Low"), 
           xlab = "Time (months)",ylab = "Disease-free survival")
STXBP1.dfs <- survfit(Surv(DFSm,DFSstate)~STXBP1.lab,
                       data = clinical.ps)
summary(STXBP1.dfs)
ggsurvplot(STXBP1.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#FF6361","#CFD3CE"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("High","Low"), 
           xlab = "Time (months)",ylab = "Disease-free survival")
CDH2.dfs <- survfit(Surv(DFSm,DFSstate)~CDH2.lab,
                      data = clinical.ps)
summary(CDH2.dfs)
ggsurvplot(CDH2.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#FF6361","#CFD3CE"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("High","Low"), 
           xlab = "Time (months)",ylab = "Disease-free survival")






















