###  Figure 4A  ###
library(ggnewscale)
library(tidyverse)
library(ggplot2)
library(ggrepel)
lcnec_ls.result <- read.csv("lcnec_vs_ls.csv",stringsAsFactors = F,
                            row.names = 1,check.names = F)
lcnec_sclc.result <- read.csv("lcnec_vs_sclc.csv",stringsAsFactors = F,
                              row.names = 1,check.names = F)
ls_sclc.result <- read.csv("ls_vs_sclc.csv",stringsAsFactors = F,
                           row.names = 1,check.names = F)
### top 10 sig
lcnec_ls.result.new <- lcnec_ls.result[c(which(lcnec_ls.result$sig == "LCNEC-up"),
                                         which(lcnec_ls.result$sig == "LS-up")),]
lcnec_ls.result.new <- lcnec_ls.result.new[order(lcnec_ls.result.new$log2FC),]
top10.sig.lcnec.ls <- lcnec_ls.result.new[c(1:5,63:67),]
top10.sig.lcnec.ls$label <- rep("lcnec vs ls",10)
lcnec_sclc.result.new <- lcnec_sclc.result[c(which(lcnec_sclc.result$sig == "LCNEC-up"),
                                             which(lcnec_sclc.result$sig == "SCLC-up")),]
lcnec_sclc.result.new <- lcnec_sclc.result.new[order(lcnec_sclc.result.new$log2FC),]
top10.sig.lcnec.sclc <- lcnec_sclc.result.new[c(1:5,135:139),]
top10.sig.lcnec.sclc$label <- rep("lcnec vs sclc",10)
ls_sclc.result.new <- ls_sclc.result[c(which(ls_sclc.result$sig == "LS-up"),
                                       which(ls_sclc.result$sig == "SCLC-up")),]
ls_sclc.result.new <- ls_sclc.result.new[order(ls_sclc.result.new$log2FC),]
top10.sig.ls.sclc <- ls_sclc.result.new[c(1:5,24:28),]
top10.sig.ls.sclc$label <- rep("ls vs sclc",10)
lcnec_ls.result.plot <- lcnec_ls.result[which(lcnec_ls.result$sig != "Non-sig"),] 
lcnec_ls.result.plot$label <- rep("lcnec vs ls",dim(lcnec_ls.result.plot)[1])
lcnec_ls.result.plot$ptname <- rownames(lcnec_ls.result.plot)
lcnec_sclc.result.plot <- lcnec_sclc.result[which(lcnec_sclc.result$sig != "Non-sig"),] 
lcnec_sclc.result.plot$label <- rep("lcnec vs sclc",dim(lcnec_sclc.result.plot)[1])
lcnec_sclc.result.plot$ptname <- rownames(lcnec_sclc.result.plot)
ls_sclc.result.plot <- ls_sclc.result[which(ls_sclc.result$sig != "Non-sig"),] 
ls_sclc.result.plot$label <- rep("ls vs sclc",dim(ls_sclc.result.plot)[1])
ls_sclc.result.plot$ptname <- rownames(ls_sclc.result.plot)

merge.plot <- rbind(lcnec_ls.result.plot,
                    lcnec_sclc.result.plot,
                    ls_sclc.result.plot)
merge.plot.top10 <- rbind(top10.sig.lcnec.ls,
                          top10.sig.lcnec.sclc,
                          top10.sig.ls.sclc)
merge.plot.top10$gene <- rownames(merge.plot.top10)

dbar <- 
  merge.plot.top10 %>% 
  group_by(label) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  select(label, log2FC_min, log2FC_max, label) %>% 
  rename(label = label)

ggplot()+
  geom_col(data = dbar, 
           mapping = aes(x = label,y = log2FC_min),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, 
           mapping = aes(x = label,y = log2FC_max),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_jitter(data=merge.plot,
              aes(x=label,y=log2FC,
                  color=sig),width =0.35) +
  geom_tile(data = merge.plot.top10, 
            aes(x = label,
                y = 0,
                fill = label),
            height=0.5,
            color = "black",
            alpha = 0.6,
            show.legend = F) +
  ggsci::scale_fill_npg() + 
  ggsci::scale_color_npg() + 
  scale_color_manual(values = c('#7A4D7B','#526187',
                                '#ED8A3F','#9B8281')) +
  ylab("log2 (Fold-change)")


###  Figure 4C left   ###
library(pheatmap)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
data.origin <- read.csv("6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
lcnec_vs_ls <- read.csv("lcnec_vs_ls.csv",stringsAsFactors = F,
                        row.names = 1,check.names = F)
lcnec_vs_sclc <- read.csv("lcnec_vs_sclc.csv",stringsAsFactors = F,
                          row.names = 1,check.names = F)
ls_vs_sclc <- read.csv("ls_vs_sclc.csv",stringsAsFactors = F,
                       row.names = 1,check.names = F)
lcnec_up <- intersect(rownames(lcnec_vs_ls[which(lcnec_vs_ls$sig == "LCNEC-up"),]),
                      rownames(lcnec_vs_sclc[which(lcnec_vs_sclc$sig == "LCNEC-up"),]))
ls_up <- intersect(rownames(lcnec_vs_ls[which(lcnec_vs_ls$sig == "LS-up"),]),
                   rownames(ls_vs_sclc[which(ls_vs_sclc$sig == "LS-up"),]))
sclc_up <- intersect(rownames(lcnec_vs_sclc[which(lcnec_vs_sclc$sig == "SCLC-up"),]),
                     rownames(ls_vs_sclc[which(ls_vs_sclc$sig == "SCLC-up"),]))
lcnec_up.all <- unique(c(rownames(lcnec_vs_ls[which(lcnec_vs_ls$sig == "LCNEC-up"),]),
                         rownames(lcnec_vs_sclc[which(lcnec_vs_sclc$sig == "LCNEC-up"),])))
ls_up.all <- unique(c(rownames(lcnec_vs_ls[which(lcnec_vs_ls$sig == "LS-up"),]),
                      rownames(ls_vs_sclc[which(ls_vs_sclc$sig == "LS-up"),])))
sclc_up.all <- unique(c(rownames(lcnec_vs_sclc[which(lcnec_vs_sclc$sig == "SCLC-up"),]),
                        rownames(ls_vs_sclc[which(ls_vs_sclc$sig == "SCLC-up"),])))
lcnec_up <- unique(lcnec_up)
ls_up <- unique(ls_up)
sclc_up <- unique(sclc_up)
intersect(lcnec_up,lcnec_up.all)
lcnec_up.delet <- unique(c(intersect(lcnec_up,ls_up.all),intersect(lcnec_up,sclc_up.all)))
ls_up.delet <- unique(c(intersect(ls_up,lcnec_up.all),intersect(ls_up,sclc_up.all)))
sclc_up.delet <- unique(c(intersect(sclc_up,lcnec_up.all),intersect(sclc_up,ls_up.all)))
lcnec_up <- lcnec_up[-match(lcnec_up.delet,lcnec_up)]
sclc_up <- sclc_up[-match(sclc_up.delet,sclc_up)]
lcnec_up.gene <- data.origin[lcnec_up,2]
ls_up.gene <- data.origin[ls_up,2]
sclc_up.gene <- data.origin[sclc_up,2]
data.origin.pro <- data.origin[c(lcnec_up,ls_up,sclc_up),]
gene <- data.origin.pro$`Gene name`
data.origin.pro <- data.origin.pro[,-c(1,2)]
data.origin.pro <- as.matrix(data.origin.pro)
rownames(data.origin.pro) <- gene
clini <- clini[colnames(data.origin.pro),]
red <- "#F44336";
blue <- "#587AA7";
yellow <- "#FFF200";
black <- "#1F1F1F";
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- apply(data.origin.pro,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data.origin.pro)
rownames(linshi) <- rownames(data.origin.pro)
hist(linshi)
linshi[linshi>1] <- 1
linshi[linshi<(-1)] <- c(-1)
annotation_col <- data.frame(type = clini$type,
                             DFSstate = as.factor(clini$DFSstate),
                             Osstate = as.factor(clini$Osstate),
                             age_60 = as.factor(clini$age_60),
                             sex = as.factor(clini$sex),
                             smoking = as.factor(clini$smoking),
                             Tumorsite = as.factor(clini$Tumor_site),
                             ajcc = as.factor(clini$AJCCstage),
                             T_stage = as.factor(clini$T_stage),
                             N_metastasis = as.factor(clini$N_metastasis),
                             M_metastasis = as.factor(clini$M_metastasis),
                             Postoperative_adjuvant_chemotherapy = as.factor(clini$Postoperative_adjuvant_chemotherapy),
                             Postoperative_adjuvant_chemotherapy_type = as.factor(clini$Postoperative_adjuvant_chemotherapy_type))
rownames(annotation_col) <- rownames(clini)
ann_colors = list(type = c("LCNEC" = "#526187","LCNEC+SCLC" = "#ED8A3F",
                           "SCLC" = "#9B8281"),
                  DFSstate = c("0" = "#A1DFDB", "1" = "#F79990"),
                  Osstate = c("0" = "#A1DFDB", "1" = "#F79990"),
                  age_60 = c("<60" = "#A1DFDB",">=60" = "#F79990"),
                  sex = c("female" = "#A1DFDB","male" = "#F79990"),
                  smoking = c("No"="#A1DFDB","Yes"="#F79990"),
                  Tumorsite = c("left" = "#8A6EAF","left-down" = "#3A4A7D","left-up" = "#7F8FA6",
                                "right" = "#FFB6C1","right-down" = "#FFCF48","right-middle" = "#a7f2a7",
                                "right-up" = "#a7a7f2"),
                  ajcc = c("I" = "#F7F2E9","II" = "#D2D2D2",
                           "III" = "#A9A9A9","IV" = "#848484"),
                  T_stage = c("1"="#dee2d1","2"="#91c0c0",
                              "3"="#c29e2e","4"="#647370"),
                  N_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  M_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy = c("No"="#A1DFDB","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy_type = c("No"="#A1DFDB","NSCLC" = "#ffad60",
                                                               "SCLC" = "#005792","SCLC+NSCLC" = "#E8222D",
                                                               "Unknown" = "#a39e9e"))
linshi <- cbind(linshi[,17:46],linshi[,1:16],linshi[,47:76])
out <- pheatmap(linshi,fontsize=6,gaps_row = c(61,70),
                color  = colorRampPalette(c(blue,white,red))(100),
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                clustering_method = "ward.D",
                border_color = "grey60",
                cluster_cols = F, cluster_rows = F,
                show_rownames = F, show_colnames = T)


###  Figure 4C  right   ###
library(ggplot2)
library(ggpubr)
enrich <- read.csv("DIP_enrichment_result_merge.csv",
                   stringsAsFactors = F,row.names = 1,check.names = F)
enrich$name <- factor(enrich$name,levels = c(enrich$name))
enrich.pl <- enrich[which(enrich$cancer_type == "LCNEC"),]
enrich.ls <- enrich[which(enrich$cancer_type == "LS"),]
enrich.ps <- enrich[which(enrich$cancer_type == "SCLC"),]

ggplot(enrich.pl, aes(x=name, y=count,fill = -log10(pvalue))) + 
  geom_bar(stat="identity",  position=position_dodge(),width=0.5) + 
  theme_bw() +
  scale_fill_gradient2(low = "#f7f7f7", high = "#8A9EB5",
                       mid = "#f7f7f7",midpoint = -log10(0.05)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=8,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black", 
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("count")
ggplot(enrich.ls, aes(x=name, y=count,fill = -log10(pvalue))) + 
  geom_bar(stat="identity",  position=position_dodge(),width=0.5) + 
  theme_bw() +
  scale_fill_gradient2(low = "#f7f7f7", high = "#ED8A3F",
                       mid = "#f7f7f7",midpoint = -log10(0.05)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=8,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("count")
ggplot(enrich.ps, aes(x=name, y=count,fill = -log10(pvalue))) + 
  geom_bar(stat="identity",  position=position_dodge(),width=0.5) + 
  theme_bw() +
  scale_fill_gradient2(low = "#f7f7f7", high = "#9B8281",
                       mid = "#f7f7f7",midpoint = -log10(0.05)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=8,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("count")





























