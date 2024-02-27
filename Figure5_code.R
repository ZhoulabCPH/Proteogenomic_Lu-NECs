
###   Figure 5A  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
map <- read.csv("ESTIMATE_result.csv",stringsAsFactors = F,
                row.names = 1,check.names = F)
map$cancer.type <- substr(rownames(map),1,2)
map$cancer.type <- factor(map$cancer.type,
                          levels = c("PL","LS","PS"))
ggplot() +
  geom_half_boxplot(data = map, 
                    aes(x=cancer.type, y=ImmuneScore, 
                        fill = cancer.type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = map, 
                  aes(x=cancer.type, y=ImmuneScore, color=cancer.type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  scale_color_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  stat_compare_means(data = map, 
                     aes(x=cancer.type, y=ImmuneScore, fill=cancer.type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="ImmuneScore")
ggplot() +
  geom_half_boxplot(data = map, 
                    aes(x=cancer.type, y=StromalScore, 
                        fill = cancer.type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = map, 
                  aes(x=cancer.type, y=StromalScore, color=cancer.type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  scale_color_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  stat_compare_means(data = map, 
                     aes(x=cancer.type, y=StromalScore, fill=cancer.type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="StromalScore")
ggplot() +
  geom_half_boxplot(data = map, 
                    aes(x=cancer.type, y=ESTIMATEScore, 
                        fill = cancer.type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = map, 
                  aes(x=cancer.type, y=ESTIMATEScore, color=cancer.type), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  scale_color_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  stat_compare_means(data = map, 
                     aes(x=cancer.type, y=ESTIMATEScore, fill=cancer.type),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(x='', y="ESTIMATEScore")


###  Figure 5B left  ###
library(pheatmap)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
imm <- read.csv("imm_cell.csv",row.names = 1,
                stringsAsFactors = F,check.names = F)
clini <- clini[colnames(imm),]
annotation_col <- data.frame(Cancer_type = as.factor(substr(rownames(clini),1,2)),
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
ann_colors = list(Cancer_type = c("PL" = "#8A9EB5","LS" = "#ED8A3F","PS" = "#9B8281"),
                  DFSstate = c("0" = "#D3E2F2", "1" = "#F79990"),
                  Osstate = c("0" = "#D3E2F2", "1" = "#F79990"),
                  age_60 = c("<60" = "#D3E2F2",">=60" = "#F79990"),
                  sex = c("female" = "#D3E2F2","male" = "#F79990"),
                  smoking = c("No"="#D3E2F2","Yes"="#F79990"),
                  Tumorsite = c("left" = "#8A6EAF","left-down" = "#3A4A7D","left-up" = "#7F8FA6",
                                "right" = "#FFB6C1","right-down" = "#FFCF48","right-middle" = "#a7f2a7",
                                "right-up" = "#a7a7f2"),
                  ajcc = c("I" = "#F7F2E9","II" = "#D2D2D2",
                           "III" = "#A9A9A9","IV" = "#848484"),
                  T_stage = c("1"="#dee2d1","2"="#91c0c0",
                              "3"="#c29e2e","4"="#647370"),
                  N_metastasis = c("No"="#D3E2F2","Yes"="#F79990"),
                  M_metastasis = c("No"="#D3E2F2","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy = c("No"="#D3E2F2","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy_type = c("No"="#D3E2F2","NSCLC" = "#ffad60",
                                                               "SCLC" = "#005792","SCLC+NSCLC" = "#E8222D",
                                                               "Unknown" = "#a39e9e"))
linshi <- apply(imm, 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(imm)
linshi[linshi < (-2)] <- (-2)
linshi[linshi > 2] <- 2
red <- "#D94E48";
blue <- "#A6A6A6";
white <- rgb(255,255,255,maxColorValue = 255)
rownames(linshi) <- rownames(imm)
out <- pheatmap(linshi,fontsize=6,
                color  = colorRampPalette(c(blue,white,red))(100),
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                clustering_method = "ward.D2",
                border_color = "grey60",
                cluster_cols = F, cluster_rows = F,
                show_rownames = T, show_colnames = T)


###  Figure 5B middle  ###
library(ggridges)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gghalves)
imm <- read.csv("imm_cell.csv",row.names = 1,
                stringsAsFactors = F,check.names = F)
linshi <- apply(imm, 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(imm)
linshi[linshi < (-2)] <- (-2)
linshi[linshi > 2] <- 2
pat <- c()
for ( i in 1:76){
  pat <- c(pat,rep(colnames(imm)[i],28))
}
###  Zscore 后免疫富集得分
imm.data <- data.frame(x = as.numeric(as.matrix(linshi)),
                       y = rep(rownames(imm),76),
                       pat = pat)
imm.data$cancer.type <- substr(imm.data$pat,1,2)
imm.data$y <- factor(imm.data$y,levels = rev(c("Activated B cell","Activated CD4 T cell",
                                               "Activated CD8 T cell","Central memory CD4 T cell",
                                               "Central memory CD8 T cell","Effector memeory CD4 T cell",
                                               "Effector memeory CD8 T cell","Gamma delta T cell",
                                               "Immature  B cell","Memory B cell",
                                               "Regulatory T cell","T follicular helper cell",
                                               "Type 1 T helper cell","Type 2 T helper cell",
                                               "Type 17 T helper cell","Natural killer T cell",
                                               "Activated dendritic cell","Eosinophil",
                                               "Immature dendritic cell","Mast cell",
                                               "Natural killer cell","Neutrophil",
                                               "Plasmacytoid dendritic cell","CD56bright natural killer cell",
                                               "CD56dim natural killer cell","Macrophage",
                                               "MDSC","Monocyte")))
ggplot(imm.data, aes(x = x, y = y,fill = cancer.type)) +
  geom_density_ridges(alpha=0.3) + 
  theme_ridges() +
  scale_fill_manual(values = c('#ed8a3f','#8a9eb5','#9b8281'))+
  xlab("Immune fraction")



###  Figure 5B right  ###
library(ggplot2)
imm <- read.csv("imm_cell.csv",check.names = F,
                stringsAsFactors = F,row.names = 1)
imm <- as.data.frame(t(imm))
imm.pl <- imm[which(substr(rownames(imm),1,2) == "PL"),]
imm.ls <- imm[which(substr(rownames(imm),1,2) == "LS"),]
imm.ps <- imm[which(substr(rownames(imm),1,2) == "PS"),]
FC_pl_vs_oth <- c()
P_pl_vs_oth <- c()
type_pl_vs_oth <- c()
FC_ls_vs_oth <- c()
P_ls_vs_oth <- c()
type_ls_vs_oth <- c()
FC_ps_vs_oth <- c()
P_ps_vs_oth <- c()
type_ps_vs_oth <- c()
P_all <- c()
for (i in 1:28){
  FC_pl_vs_oth <- c(FC_pl_vs_oth,(mean(imm.pl[,i]) / mean(rbind(imm.ls,imm.ps)[,i])))
  P_pl_vs_oth <- c(P_pl_vs_oth,wilcox.test(imm.pl[,i],rbind(imm.ls,imm.ps)[,i])$p.value)
  FC_ls_vs_oth <- c(FC_ls_vs_oth,(mean(imm.ls[,i]) / mean(rbind(imm.pl,imm.ps)[,i])))
  P_ls_vs_oth <- c(P_ls_vs_oth,wilcox.test(imm.ls[,i],rbind(imm.pl,imm.ps)[,i])$p.value)
  FC_ps_vs_oth <- c(FC_ps_vs_oth,(mean(imm.ps[,i]) / mean(rbind(imm.pl,imm.ls)[,i])))
  P_ps_vs_oth <- c(P_ps_vs_oth,wilcox.test(imm.ps[,i],rbind(imm.pl,imm.ls)[,i])$p.value)
  P_all <- c(P_all,kruskal.test(imm[,i]~substr(rownames(imm),1,2))$p.value)
}
result.mat <- data.frame(fc = c(FC_pl_vs_oth,FC_ls_vs_oth,FC_ps_vs_oth),
                         p = c(P_pl_vs_oth,P_ls_vs_oth,P_ps_vs_oth),
                         type = c(rep("pl_vs_oth",28),rep("ls_vs_oth",28),
                                  rep("ps_vs_oth",28)),
                         imm = colnames(imm.pl))
result.mat$type <- factor(result.mat$type,
                          levels = c("pl_vs_oth","ls_vs_oth","ps_vs_oth"))
result.mat.pl <- result.mat[which(result.mat$type == "pl_vs_oth"),]
result.mat.ls <- result.mat[which(result.mat$type == "ls_vs_oth"),]
result.mat.ps <- result.mat[which(result.mat$type == "ps_vs_oth"),]
imm <- read.csv("imm_cell.csv",row.names = 1,
                stringsAsFactors = F,check.names = F)
linshi <- apply(imm, 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(imm)
linshi[linshi < (-2)] <- (-2)
linshi[linshi > 2] <- 2
linshi.pl <- linshi[result.mat.pl$imm,
                    which(substr(colnames(linshi),1,2) == "PL")]
linshi.ls <- linshi[result.mat.ls$imm,
                    which(substr(colnames(linshi),1,2) == "LS")]
linshi.ps <- linshi[result.mat.ps$imm,
                    which(substr(colnames(linshi),1,2) == "PS")]
p.mat <- data.frame(p_pl = result.mat.pl$p,
                    p_ls = result.mat.ls$p,
                    p_ps = result.mat.ps$p,
                    p_all = P_all)
rownames(p.mat) <- result.mat.pl$imm
p.mat.dot <- data.frame(P = c(result.mat.pl$p,result.mat.ls$p,
                              result.mat.ps$p,P_all),
                        immcell = rep(rownames(p.mat),4),
                        type = c(rep("pl",28),rep("ls",28),
                                 rep("ps",28),rep("all",28)))
p.mat.dot$type <- factor(p.mat.dot$type,
                         levels = c("all","pl","ls","ps"))
p.mat.dot$immcell <- factor(p.mat.dot$immcell,
                            levels = result.mat.pl$imm)
ggplot(data = p.mat.dot, mapping = aes(x=type,y=immcell))+
         geom_point(aes(size= -log10(P)),color= "black")
p.mat1 <- data.frame(p_pl = result.mat.pl$p,
                    p_ls = result.mat.ls$p,
                    p_ps = result.mat.ps$p,
                    inten_pl = rowMeans(linshi.pl),
                    inten_ls = rowMeans(linshi.ls),
                    inten_ps = rowMeans(linshi.ps))
rownames(p.mat1) <- result.mat.pl$imm
p.mat.dot1 <- data.frame(P = c(result.mat.pl$p,result.mat.ls$p,
                              result.mat.ps$p),
                        immcell = rep(rownames(p.mat),3),
                        intensity = c(p.mat1$inten_pl,
                                      p.mat1$inten_ls,
                                      p.mat1$inten_ps),
                        type = c(rep("pl",28),rep("ls",28),
                                 rep("ps",28)))
p.mat.dot1$type <- factor(p.mat.dot1$type,
                         levels = c("pl","ls","ps"))
p.mat.dot1$immcell <- factor(p.mat.dot1$immcell,
                            levels = result.mat.pl$imm)
ggplot(data = p.mat.dot1, mapping = aes(x=type,y=immcell))+
  geom_point(aes(size= -log10(P),color= intensity)) +
  scale_color_gradient2(low = "#A6A6A6",mid = "#F9F7F7",
                        high = "#D94E48",midpoint = 0)



##  Figure 5C  ###
hall <- read.csv("hallmark.csv",check.names = F,
                 stringsAsFactors = F,row.names = 1)
fge <- read.table("29_FGE.txt",check.names = F,
                  stringsAsFactors = F,row.names = 1,sep="\t")
fge <- fge[,colnames(hall)]
hall.t <- apply(hall, 1, scale)
hall.t <- t(hall.t)
colnames(hall.t) <- colnames(hall)
fge.t <- apply(fge, 1, scale)
fge.t <- t(fge.t)
colnames(fge.t) <- colnames(fge)
all <- rbind(hall.t,fge.t)
all[all > 2] <- 2
all[all < (-2)] <- (-2)
all <- t(all)
all <- rbind(all[17:46,],all[1:16,],all[47:76,])
yellow <- "#FFF200";
white <- "#FFFFFF";
purple <- "#844391"
pheatmap(all,fontsize=6,cellwidth = 6,cellheight = 2,
         gaps_row = c(30,46),
         color  = colorRampPalette(c(purple,white,yellow))(100),
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T
)  



























