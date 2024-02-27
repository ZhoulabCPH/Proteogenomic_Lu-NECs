###  Figure 2A  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("snvInDel_zzc.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
data_all <- data.frame(y = c(data$`Total SNVs`/42,data$InDels/42,data$MNV/42),
                       vari_type = c(rep("SNV",dim(data)[1]),rep("Indels",dim(data)[1]),rep("MNV",dim(data)[1])),
                       cancer_type = c(type$type,type$type,type$type),
                       patient = c(rownames(data),rownames(data),rownames(data)))
data_all$patient <- factor(data_all$patient,levels = rownames(data))
ggplot(data_all,mapping = aes(patient,y,fill=vari_type))+
  geom_bar(stat='identity',position='stack') +
  labs(x = 'patient',y = 'Mutations per Mb') +
  scale_fill_manual(values = c("#ffd800","#fad3cf","#a1bad0"))+
  theme(axis.title =element_text(size = 8),
        axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
data$TMB <- c((data$`Total SNVs`/42)+(data$InDels/42)+(data$MNV/42))
data$Cancer_type <- substr(data$patient,1,2)
data$Cancer_type <- factor(data$Cancer_type,
                           levels = c("PL","LS","PS"))
ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Cancer_type, y=TMB, 
                        fill = Cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Cancer_type, y=TMB, color=Cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data, 
                     aes(x=Cancer_type, y=TMB, fill=Cancer_type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="TMB")
tmb <- data.frame(y = c((data$`Total SNVs`/42)+(data$InDels/42)+(data$MNV/42)),
                  vari_type = c(rep("TMB",dim(data)[1])),
                  patient = c(rownames(data)))
#write.csv(tmb,"TMB.csv",quote = F)
data.pl <- data[which(data$Cancer_type == "PL"),]
data.ls <- data[which(data$Cancer_type == "LS"),]
data.ps <- data[which(data$Cancer_type == "PS"),]
median(data.pl$TMB)
median(data.ls$TMB)
median(data.ps$TMB)
wilcox.test(data.ls$TMB,data.ps$TMB)
#####
### clinical Figure 2A ###
library(pheatmap)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
data <- read.csv("snvInDel_zzc.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
clini <- clini[rownames(data),]
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
linshi <- t(data.frame(x = 1:93, y = 1:93))
colnames(linshi) <- rownames(clini)
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
pheatmap(linshi,fontsize=6,
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         clustering_method = "ward.D2",
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = T)


###  Figure 2B  ###
library(plyr)
library(ggplot2)
library(ggpubr)
ID <- read.csv("ID.csv",
               stringsAsFactors = F,row.names = 1,check.names = F)
maf <- matrix(,,3)
colnames(maf) <- c("sample","Variant_Classification","Variant_Type")
for (i in 1:dim(ID)[1]){
  a <- paste(rownames(ID)[i],"_TN.PASS.anno.filtered.maf",sep = "")
  data <- read.table(a)
  colnames(data) <- data[1,]
  data <- data[-1,]
  type <- table(data$Variant_Type)
  maf <- rbind(maf,
               cbind(rep(ID$patient[i],dim(data)[1]),
                     data$Variant_Classification,
                     data$Variant_Type))
}
maf <- as.data.frame(maf[-1,])
data <- read.csv("snvInDel_zzc.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
maf$sample <- factor(maf$sample,levels = c(rownames(data)))
maf <- maf[which(maf$Variant_Classification != "Translation_Start_Site"),]
maf <- maf[which(maf$Variant_Classification != "Inframe_INDEL"),]
maf$number <- 1
data_duiji1 <- ddply(maf,'sample',transform,percent = 1/sum(number)*100)
ggplot(data_duiji1)+
  scale_fill_manual(values = c("#C7DCE8","#1E78B0","#FFFDCF","#FFC857",
                               "#AEE18A","#C5B3D4","#F0BDBD",
                               "#E61C19","#FF8A00","#a39e9e"))+ 
  geom_bar(aes(x=sample,fill=Variant_Classification),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())
data_duiji1_PL <- data_duiji1[which(substr(data_duiji1$sample,1,2) == "PL"),]
data_duiji1_LS <- data_duiji1[which(substr(data_duiji1$sample,1,2) == "LS"),]
data_duiji1_PS <- data_duiji1[which(substr(data_duiji1$sample,1,2) == "PS"),]
table(data_duiji1_PL$Variant_Classification)
table(data_duiji1_LS$Variant_Classification)
table(data_duiji1_PS$Variant_Classification)
table(data_duiji1_PL$Variant_Classification)/sum(table(data_duiji1_PL$Variant_Classification))
table(data_duiji1_LS$Variant_Classification)/sum(table(data_duiji1_LS$Variant_Classification))
table(data_duiji1_PS$Variant_Classification)/sum(table(data_duiji1_PS$Variant_Classification))



###  Figure 2C ###
library(plyr)
library(ggplot2)
library(ggpubr)
ID <- read.csv("ID.csv",
               stringsAsFactors = F,row.names = 1,check.names = F)
maf <- matrix(,,4)
colnames(maf) <- c("sample","Variant_Classification","Variant_Type","Mutation_changes")
for (i in 1:dim(ID)[1]){
  a <- paste(rownames(ID)[i],"_TN.PASS.anno.filtered.maf",sep = "")
  data <- read.table(a)
  colnames(data) <- data[1,]
  data <- data[-1,]
  type <- table(data$Variant_Type)
  maf <- rbind(maf,
               cbind(rep(ID$patient[i],dim(data)[1]),
                     data$Variant_Classification,
                     data$Variant_Type,
                     paste(data$Reference_Allele,data$Tumor_Seq_Allele2,sep = ">")))
}
maf <- as.data.frame(maf[-1,])
data <- read.csv("snvInDel_zzc.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
maf$sample <- factor(maf$sample,levels = c(rownames(data)))
maf_snp <- maf[which(maf$Variant_Type == "SNP"),]
maf_snp[maf_snp$Mutation_changes == "G>T",4] <- "C>A"
maf_snp[maf_snp$Mutation_changes == "G>C",4] <- "C>G"
maf_snp[maf_snp$Mutation_changes == "G>A",4] <- "C>T"
maf_snp[maf_snp$Mutation_changes == "A>T",4] <- "T>A"
maf_snp[maf_snp$Mutation_changes == "A>G",4] <- "T>C"
maf_snp[maf_snp$Mutation_changes == "A>C",4] <- "T>G"
maf_dnp <- maf[which(maf$Variant_Type == "DNP"),]
maf_tnp <- maf[which(maf$Variant_Type == "TNP"),]
maf_snp$number <- 1
maf_snp <- ddply(maf_snp,'sample',transform,percent = 1/sum(number)*100)
ggplot(maf_snp)+
  scale_fill_manual(values = c("#A9BCD0","#000000","#E83F3D","#C7C7C7","#A1CC63","#E8C2BF"))+ 
  geom_bar(aes(x=sample,fill=Mutation_changes),position = "fill")+
  theme_bw()+ theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"),
        legend.text=element_text(colour="black",size=10),
        legend.title=element_text(colour="black", size=10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
data_new <- matrix(,93,6)
colnames(data_new) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
for (i in 1:dim(type)[1]){
  weizhi <- which(rownames(type)[i] == maf_snp$sample)
  maf_snp_new <- maf_snp[weizhi,]
  data_new[i,] <- c(length(which(maf_snp_new$Mutation_changes == "C>A")),
                    length(which(maf_snp_new$Mutation_changes == "C>G")),
                    length(which(maf_snp_new$Mutation_changes == "C>T")),
                    length(which(maf_snp_new$Mutation_changes == "T>A")),
                    length(which(maf_snp_new$Mutation_changes == "T>C")),
                    length(which(maf_snp_new$Mutation_changes == "T>G")))
}
rownames(data_new) <- rownames(type)
data_new <- as.data.frame(data_new)
data_new$Cancer_type <- type$type
data_new1 <- data.frame(y = c(data_new$`C>A`,data_new$`C>G`,data_new$`C>T`,
                              data_new$`T>A`,data_new$`T>C`,data_new$`T>G`),
                        change_type = c(rep("C>A",93),rep("C>G",93),rep("C>T",93),
                                        rep("T>A",93),rep("T>C",93),rep("T>G",93)),
                        Cancer_type = c(data_new$Cancer_type,data_new$Cancer_type,data_new$Cancer_type,
                                        data_new$Cancer_type,data_new$Cancer_type,data_new$Cancer_type))
ggplot() +
  geom_half_boxplot(data = data_new1, 
                    aes(x=change_type, y=y, 
                        fill = Cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data_new1, 
                  aes(x=change_type, y=y, color=Cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data_new1, 
                     aes(x=change_type, y=y, fill=Cancer_type),
                     method = "kruskal.test")+   theme_bw()+
  scale_y_continuous(expand = c(0,0), trans = scales::pseudo_log_trans(), 
                     breaks = c(5,10, 30, 100, 250, 600, 1500,3000),
                     limits = c(5, 3000))+
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             alpha=1) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")
data_new_l <- data_new[which(data_new$Cancer_type == "LCNEC"),]
data_new_ls <- data_new[which(data_new$Cancer_type == "LCNEC+SCLC"),]
data_new_s <- data_new[which(data_new$Cancer_type == "SCLC"),]
wilcox.test(data_new_l$`C>A`,data_new_ls$`C>A`)
wilcox.test(data_new_l$`C>A`,data_new_s$`C>A`)
wilcox.test(data_new_ls$`C>A`,data_new_s$`C>A`)
median(data_new_l$`T>G`)
median(data_new_ls$`T>G`)
median(data_new_s$`T>G`)
data.ti.tv <- data.frame(titv = (data_new$`C>T`+
                                   data_new$`T>C`)/(data_new$`C>A`+
                                                      data_new$`C>G`+
                                                      data_new$`T>A`+
                                                      data_new$`T>G`),
                         Cancer_type = data_new$Cancer_type,
                         patient = rownames(data_new))
ggplot() +
  geom_half_boxplot(data = data.ti.tv, 
                    aes(x=Cancer_type, y=titv, 
                        fill = Cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data.ti.tv, 
                  aes(x=Cancer_type, y=titv, color=Cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data.ti.tv, 
                     aes(x=Cancer_type, y=titv, fill=Cancer_type),
                     method = "kruskal.test")+   theme_bw()+
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             alpha=1) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")
data.ti.tv_l <- data.ti.tv[which(data.ti.tv$Cancer_type == "LCNEC"),]
data.ti.tv_ls <- data.ti.tv[which(data.ti.tv$Cancer_type == "LCNEC+SCLC"),]
data.ti.tv_s <- data.ti.tv[which(data.ti.tv$Cancer_type == "SCLC"),]
wilcox.test(data.ti.tv_l$titv,data.ti.tv_ls$titv)
wilcox.test(data.ti.tv_l$titv,data.ti.tv_s$titv)
wilcox.test(data.ti.tv_ls$titv,data.ti.tv_s$titv)
median(data.ti.tv_l$titv)
median(data.ti.tv_ls$titv)
median(data.ti.tv_s$titv)


###  Figure 2D ###
library(plyr)
library(ggplot2)
library(ggpubr)
ID <- read.csv("ID.csv",
               stringsAsFactors = F,row.names = 1,check.names = F)
maf <- matrix(,,4)
colnames(maf) <- c("sample","Variant_Classification","Variant_Type","Mutation_changes")
for (i in 1:dim(ID)[1]){
  a <- paste(rownames(ID)[i],"_TN.PASS.anno.filtered.maf",sep = "")
  data <- read.table(a)
  colnames(data) <- data[1,]
  data <- data[-1,]
  type <- table(data$Variant_Type)
  maf <- rbind(maf,
               cbind(rep(ID$patient[i],dim(data)[1]),
                     data$Variant_Classification,
                     data$Variant_Type,
                     paste(data$Reference_Allele,data$Tumor_Seq_Allele2,sep = ">")))
}
maf <- as.data.frame(maf[-1,])
data <- read.csv("snvInDel_zzc.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
maf$sample <- factor(maf$sample,levels = c(rownames(data)))
maf_dnp <- maf[which(maf$Variant_Type == "DNP"),]

maf_dnp[c(which(maf_dnp$Mutation_changes == "AC>CA"),which(maf_dnp$Mutation_changes == "AC>CG"),
          which(maf_dnp$Mutation_changes == "AC>CT"),which(maf_dnp$Mutation_changes == "AC>GA"),
          which(maf_dnp$Mutation_changes == "AC>GG"),which(maf_dnp$Mutation_changes == "AC>GT"),
          which(maf_dnp$Mutation_changes == "AC>TA"),which(maf_dnp$Mutation_changes == "AC>TG"),
          which(maf_dnp$Mutation_changes == "AC>TT")),4] <- "AC>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "AT>CA"),which(maf_dnp$Mutation_changes == "AT>CC"),
          which(maf_dnp$Mutation_changes == "AT>CG"),which(maf_dnp$Mutation_changes == "AT>GA"),
          which(maf_dnp$Mutation_changes == "AT>GC"),
          which(maf_dnp$Mutation_changes == "AT>TA")),4] <- "AT>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "CC>AA"),which(maf_dnp$Mutation_changes == "CC>AG"),
          which(maf_dnp$Mutation_changes == "CC>AT"),which(maf_dnp$Mutation_changes == "CC>GA"),
          which(maf_dnp$Mutation_changes == "CC>GG"),which(maf_dnp$Mutation_changes == "CC>GT"),
          which(maf_dnp$Mutation_changes == "CC>TA"),which(maf_dnp$Mutation_changes == "CC>TG"),
          which(maf_dnp$Mutation_changes == "CC>TT")),4] <- "CC>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "CC>TT")),4] <- "CC>TT"
maf_dnp[c(which(maf_dnp$Mutation_changes == "CG>AT"),which(maf_dnp$Mutation_changes == "CG>GC"),
          which(maf_dnp$Mutation_changes == "CG>GT"),which(maf_dnp$Mutation_changes == "CG>TA"),
          which(maf_dnp$Mutation_changes == "CG>TC"),
          which(maf_dnp$Mutation_changes == "CG>TT")),4] <- "CG>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "CT>AA"),which(maf_dnp$Mutation_changes == "CT>AC"),
          which(maf_dnp$Mutation_changes == "CT>AG"),which(maf_dnp$Mutation_changes == "CT>GA"),
          which(maf_dnp$Mutation_changes == "CT>GC"),which(maf_dnp$Mutation_changes == "CT>GG"),
          which(maf_dnp$Mutation_changes == "CT>TA"),which(maf_dnp$Mutation_changes == "CT>TC"),
          which(maf_dnp$Mutation_changes == "CT>TG")),4] <- "CT>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "GC>AA"),which(maf_dnp$Mutation_changes == "GC>AG"),
          which(maf_dnp$Mutation_changes == "GC>AT"),which(maf_dnp$Mutation_changes == "GC>CA"),
          which(maf_dnp$Mutation_changes == "GC>CG"),
          which(maf_dnp$Mutation_changes == "GC>TA")),4] <- "GC>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "TA>AT"),which(maf_dnp$Mutation_changes == "TA>CG"),
          which(maf_dnp$Mutation_changes == "TA>CT"),which(maf_dnp$Mutation_changes == "TA>GC"),
          which(maf_dnp$Mutation_changes == "TA>GG"),
          which(maf_dnp$Mutation_changes == "TA>GT")),4] <- "TA>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "TC>AA"),which(maf_dnp$Mutation_changes == "TC>AG"),
          which(maf_dnp$Mutation_changes == "TC>AT"),which(maf_dnp$Mutation_changes == "TC>CA"),
          which(maf_dnp$Mutation_changes == "TC>CG"),which(maf_dnp$Mutation_changes == "TC>CT"),
          which(maf_dnp$Mutation_changes == "TC>GA"),which(maf_dnp$Mutation_changes == "TC>GG"),
          which(maf_dnp$Mutation_changes == "TC>GT")),4] <- "TC>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "TG>AA"),which(maf_dnp$Mutation_changes == "TG>AC"),
          which(maf_dnp$Mutation_changes == "TG>AT"),which(maf_dnp$Mutation_changes == "TG>CA"),
          which(maf_dnp$Mutation_changes == "TG>CC"),which(maf_dnp$Mutation_changes == "TG>CT"),
          which(maf_dnp$Mutation_changes == "TG>GA"),which(maf_dnp$Mutation_changes == "TG>GC"),
          which(maf_dnp$Mutation_changes == "TG>GT")),4] <- "TG>NN"
maf_dnp[c(which(maf_dnp$Mutation_changes == "TT>AA"),which(maf_dnp$Mutation_changes == "TT>AC"),
          which(maf_dnp$Mutation_changes == "TT>AG"),which(maf_dnp$Mutation_changes == "TT>CA"),
          which(maf_dnp$Mutation_changes == "TT>CC"),which(maf_dnp$Mutation_changes == "TT>CG"),
          which(maf_dnp$Mutation_changes == "TT>GA"),which(maf_dnp$Mutation_changes == "TT>GC"),
          which(maf_dnp$Mutation_changes == "TT>GG")),4] <- "TT>NN"
maf_dnp <- maf_dnp[c(which(maf_dnp$Mutation_changes == "AC>NN"),
                     which(maf_dnp$Mutation_changes == "AT>NN"),
                     which(maf_dnp$Mutation_changes == "CC>NN"),
                     which(maf_dnp$Mutation_changes == "CC>TT"),
                     which(maf_dnp$Mutation_changes == "CG>NN"),
                     which(maf_dnp$Mutation_changes == "CT>NN"),
                     which(maf_dnp$Mutation_changes == "GC>NN"),
                     which(maf_dnp$Mutation_changes == "TA>NN"),
                     which(maf_dnp$Mutation_changes == "TC>NN"),
                     which(maf_dnp$Mutation_changes == "TG>NN"),
                     which(maf_dnp$Mutation_changes == "TT>NN")),]
maf_dnp$number <- 1
maf_dnp <- ddply(maf_dnp,'sample',transform,percent = 1/sum(number)*100)
ggplot(maf_dnp)+
  scale_fill_manual(values = c("#C7DCE8","#1E78B0","#85CCEB","#AEE18A","#FFFDCF","#F0BDBD","#E61C19","#FDBA6E","#F07D17","#C5B3D4"))+ #设置填充的颜色
  geom_bar(aes(x=sample,fill=Mutation_changes),position = "fill")+
  theme_bw()+ theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 45), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"),
        legend.text=element_text(colour="black",size=10),
        legend.title=element_text(colour="black", size=10),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
data_new <- matrix(,93,11)
colnames(data_new) <- c("AC>NN","AT>NN","CC>NN","CC>TT","CG>NN","CT>NN",
                        "GC>NN","TA>NN","TC>NN","TG>NN","TT>NN")
for (i in 1:dim(type)[1]){
  weizhi <- which(rownames(type)[i] == maf_dnp$sample)
  maf_snp_new <- maf_dnp[weizhi,]
  data_new[i,] <- c(length(which(maf_snp_new$Mutation_changes == "AC>NN")),
                    length(which(maf_snp_new$Mutation_changes == "AT>NN")),
                    length(which(maf_snp_new$Mutation_changes == "CC>NN")),
                    length(which(maf_snp_new$Mutation_changes == "CC>TT")),
                    length(which(maf_snp_new$Mutation_changes == "CG>NN")),
                    length(which(maf_snp_new$Mutation_changes == "CT>NN")),
                    length(which(maf_snp_new$Mutation_changes == "GC>NN")),
                    length(which(maf_snp_new$Mutation_changes == "TA>NN")),
                    length(which(maf_snp_new$Mutation_changes == "TC>NN")),
                    length(which(maf_snp_new$Mutation_changes == "TG>NN")),
                    length(which(maf_snp_new$Mutation_changes == "TT>NN")))
}
rownames(data_new) <- rownames(type)
data_new <- as.data.frame(data_new)
data_new$Cancer_type <- type$type
data_new1 <- data.frame(y = c(data_new$`AC>NN`,data_new$`AT>NN`,data_new$`CC>NN`,
                              data_new$`CC>TT`,data_new$`CG>NN`,data_new$`CT>NN`,
                              data_new$`GC>NN`,data_new$`TA>NN`,data_new$`TC>NN`,
                              data_new$`TG>NN`,data_new$`TT>NN`),
                        change_type = c(rep("AC>NN",93),rep("AT>NN",93),rep("CC>NN",93),
                                        rep("CC>TT",93),rep("CG>NN",93),rep("CT>NN",93),
                                        rep("GC>NN",93),rep("TA>NN",93),rep("TC>NN",93),
                                        rep("TG>NN",93),rep("TT>NN",93)),
                        Cancer_type = c(data_new$Cancer_type,data_new$Cancer_type,data_new$Cancer_type,
                                        data_new$Cancer_type,data_new$Cancer_type,data_new$Cancer_type,
                                        data_new$Cancer_type,data_new$Cancer_type,data_new$Cancer_type,
                                        data_new$Cancer_type,data_new$Cancer_type))

ggplot() +
  geom_half_boxplot(data = data_new1, 
                    aes(x=change_type, y=y, 
                        fill = Cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data_new1, 
                  aes(x=change_type, y=y, color=Cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data_new1, 
                     aes(x=change_type, y=y),
                     method = "kruskal.test")+   theme_bw()+
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             alpha=1) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.6),
        legend.position="right")
data_new_l <- data_new[which(data_new$Cancer_type == "LCNEC"),]
data_new_ls <- data_new[which(data_new$Cancer_type == "LCNEC+SCLC"),]
data_new_s <- data_new[which(data_new$Cancer_type == "SCLC"),]
median(data_new_l$`TT>NN`)
median(data_new_ls$`TT>NN`)
median(data_new_s$`TT>NN`)


### Figure 2E  ####
library(ggplot2)
library(ggpubr)
library(gghalves)
laml.sig.contribution <- read.csv("cosmic.csv",
                                  stringsAsFactors = F,row.names = 1,check.names = F)
laml.sig.contribution_data <- data.frame(patient = rep(colnames(laml.sig.contribution),6),
                                         y = c(as.numeric(laml.sig.contribution[1,]),
                                               as.numeric(laml.sig.contribution[2,]),
                                               as.numeric(laml.sig.contribution[3,]),
                                               as.numeric(laml.sig.contribution[4,]),
                                               as.numeric(laml.sig.contribution[5,]),
                                               as.numeric(laml.sig.contribution[6,])),
                                         cosmic = c(rep(rownames(laml.sig.contribution)[1],93),
                                                    rep(rownames(laml.sig.contribution)[2],93),
                                                    rep(rownames(laml.sig.contribution)[3],93),
                                                    rep(rownames(laml.sig.contribution)[4],93),
                                                    rep(rownames(laml.sig.contribution)[5],93),
                                                    rep(rownames(laml.sig.contribution)[6],93)))
data <- read.csv("snvInDel_zzc.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
laml.sig.contribution_data$patient <- factor(laml.sig.contribution_data$patient,
                                             levels = c(rownames(data)))
ggplot(laml.sig.contribution_data,mapping = aes(patient,y,fill=cosmic))+
  geom_bar(stat='identity',position='stack') +
  labs(x = 'patient',y = 'Cosmic signure persentage') +
  scale_fill_manual(values = c("#C69E87","#D9C0AD","#B9B7A4",
                               "#8CB09E","#6C705E","#A4A4A4"))+
  theme(axis.title =element_text(size = 8),
        axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

data_new1 <- laml.sig.contribution_data
data_new1$cancertype <- substr(data_new1$patient,1,2)
data_new1$cancertype <- factor(data_new1$cancertype,
                               levels = c("PL","LS","PS"))
data_new1$y1 <- round(data_new1$y,1)
data_new1 <- data_new1[which(data_new1$cosmic != "Unknown aetiology"),]
ggplot() +
  geom_half_boxplot(data = data_new1, 
                    aes(x=cosmic, y=y1, 
                        fill = cancertype),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  ylim(0,1) +
  geom_half_point(data = data_new1, 
                  aes(x=cosmic, y=y1, color=cancertype), 
                  shape = 16, alpha=1, size = 0.8,
                  transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#8A9EB5','#F1A66D','#B9A8A7'))+
  scale_color_manual(values = c('#8A9EB5','#F1A66D','#B9A8A7'))+
  stat_compare_means(data = data_new1, 
                     aes(x=cosmic, y=y1),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="Cosmic SBS signature (v3)")
sbs13 <- data_new1[which(data_new1$cosmic == "APOBEC Cytidine Deaminase (SBS.13)"),]
sbs3 <- data_new1[which(data_new1$cosmic == "Defects in DNA-DSB repair by HR (SBS.3)"),]
sbs6 <- data_new1[which(data_new1$cosmic == "Defective DNA mismatch repair (SBS.6)"),]
sbs4 <- data_new1[which(data_new1$cosmic == "Exposure to tobacco (smoking) mutagens (SBS.4)"),]
sbs25 <- data_new1[which(data_new1$cosmic == "Possible exposure to chemotherapy (SBS.25)"),]
median(sbs25[which(sbs25$cancertype == "PL"),2])
median(sbs25[which(sbs25$cancertype == "LS"),2])
median(sbs25[which(sbs25$cancertype == "PS"),2])


###  Figure 2F  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("snvInDel_zzc.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
wgii <- read.csv("wGII.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
wgii$patient <- rownames(wgii)
wgii_pl <- wgii[rownames(pl),]
wgii_ls <- wgii[rownames(ls),]
wgii_ps <- wgii[rownames(ps),]
data1 <- as.data.frame(rbind(wgii_pl,wgii_ls,wgii_ps))
data_all <- data.frame(y = c(wgii$wGII),
                       vari_type = c(rep("wgii",dim(data1)[1])),
                       cancer_type = c(type$type),
                       patient = c(rownames(data1)))
data_all$patient <- factor(data_all$patient,levels = rownames(data))
ggplot(data_all,aes(patient,y))+
  geom_bar(aes(fill=y),stat='identity') +
  labs(x = 'patient',y = 'wGII') +   
  scale_fill_gradient(low = "#EFEFEF", high = "#FFC2B3", na.value = NA) +
  theme(axis.title =element_text(size = 8),
        axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot() +
  geom_half_boxplot(data = data_all, 
                    aes(x=cancer_type, y=y, 
                        fill = cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data_all, 
                  aes(x=cancer_type, y=y, color=cancer_type), 
                  shape = 16, alpha=1, size = 0.8,
                  transformation = position_jitter(height = 0))+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data_all, 
                     aes(x=cancer_type, y=y),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="wGII")
median(data_all[which(data_all$cancer_type == "LCNEC"),1])
median(data_all[which(data_all$cancer_type == "LCNEC+SCLC"),1])
median(data_all[which(data_all$cancer_type == "SCLC"),1])
wilcox.test(data_all[which(data_all$cancer_type == "LCNEC+SCLC"),1],
            data_all[which(data_all$cancer_type == "SCLC"),1])


###  Figure 2G  ###
library(plyr)
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("snvInDel_zzc.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
data1 <- read.csv("MSIscore.csv",
                  stringsAsFactors = F,row.names = 1,check.names = F)
data1$patient <- rownames(data1)
data_pl <- data1[rownames(pl),]
data_ls <- data1[rownames(ls),]
data_ps <- data1[rownames(ps),]
data1 <- as.data.frame(rbind(data_pl,data_ls,data_ps))
data_all <- data.frame(y = c(data1$`MSIscore%`),
                       vari_type = c(rep("MSIscore%",dim(data1)[1])),
                       cancer_type = c(type$type),
                       patient = c(rownames(data1)))
data_all$patient <- factor(data_all$patient,levels = rownames(data))
ggplot(data_all,aes(patient,y))+
  geom_bar(aes(fill=y),stat='identity') +
  labs(x = 'patient',y = 'MSI score %') +  ylim(0,4) + 
  scale_fill_gradient(low = "#EFEFEF", high = "#C6B497", na.value = NA) +
  theme(axis.title =element_text(size = 8),
        axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
data_all$cancer_type  <- factor(data_all$cancer_type,levels = c("LCNEC",'LCNEC+SCLC','SCLC'))
ggplot() +
  geom_half_boxplot(data = data_all, 
                    aes(x=cancer_type, y=y, 
                        fill = cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data_all, 
                  aes(x=cancer_type, y=y, color=cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data_all, 
                     aes(x=cancer_type, y=y, fill=cancer_type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="MSIscore %")
median(data_all[which(data_all$cancer_type == "LCNEC"),1])
median(data_all[which(data_all$cancer_type == "LCNEC+SCLC"),1])
median(data_all[which(data_all$cancer_type == "SCLC"),1])
wilcox.test(data_all[which(data_all$cancer_type == "LCNEC+SCLC"),1],
            data_all[which(data_all$cancer_type == "SCLC"),1])

###  Figure 2H  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("snvInDel_zzc.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
data$patient <- rownames(data)
type <- read.csv("type.csv",
                 stringsAsFactors = F,
                 row.names = 1,check.names = F)
pl <- type[which(type$type == "LCNEC"),]
ls <- type[which(type$type == "LCNEC+SCLC"),]
ps <- type[which(type$type == "SCLC"),]
data_pl <- data[rownames(pl),]
data_ls <- data[rownames(ls),]
data_ps <- data[rownames(ps),]
data_pl <- data_pl[rev(order(data_pl$`Total SNVs/InDels`)),]
data_ls <- data_ls[rev(order(data_ls$`Total SNVs/InDels`)),]
data_ps <- data_ps[rev(order(data_ps$`Total SNVs/InDels`)),]
data <- as.data.frame(rbind(data_pl,data_ls,data_ps))
cnv.burden <- read.csv("CNV_burden.csv",stringsAsFactors = F,
                       row.names = 1,check.names = F)
cnv.burden$patient <- rownames(cnv.burden)
cnv.burden_pl <- cnv.burden[rownames(pl),]
cnv.burden_ls <- cnv.burden[rownames(ls),]
cnv.burden_ps <- cnv.burden[rownames(ps),]
data1 <- as.data.frame(rbind(cnv.burden_pl,cnv.burden_ls,cnv.burden_ps))
data_all <- data.frame(y = c(cnv.burden$CNV_burden),
                       vari_type = c(rep("cnv.burden",dim(data1)[1])),
                       cancer_type = c(type$type),
                       patient = c(rownames(data1)))
data_all$patient <- factor(data_all$patient,levels = rownames(data))
ggplot(data_all,aes(patient,y))+
  geom_bar(aes(fill=y),stat='identity') +
  labs(x = 'patient',y = 'CNV burden') +  ylim(0,1) + 
  scale_fill_gradient(low = "#EFEFEF", high = "#B2D7D2", na.value = NA) +
  theme(axis.title =element_text(size = 8),
        axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot() +
  geom_half_boxplot(data = data_all, 
                    aes(x=cancer_type, y=y, 
                        fill = cancer_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data_all, 
                  aes(x=cancer_type, y=y, color=cancer_type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  scale_color_manual(values = c('#526187','#ED8A3F','#9B8281'))+
  stat_compare_means(data = data_all, 
                     aes(x=cancer_type, y=y, fill=cancer_type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="CNV burden")
median(data_all[which(data_all$cancer_type == "LCNEC"),1])
median(data_all[which(data_all$cancer_type == "LCNEC+SCLC"),1])
median(data_all[which(data_all$cancer_type == "SCLC"),1])

data_all.new <- data_all[c(which(data_all$cancer_type == "LCNEC"),
                           which(data_all$cancer_type == "SCLC")),]
wilcox.test(data_all.new$y~data_all.new$cancer_type)















