###  Figure 1B left  ###
library(circlize)
library(ggplot2)
library(ggforce)
clinical <- read.csv("type.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- as.matrix(clinical)
split <- factor(substr(rownames(clinical),1,2),
                levels = c("PL","LS","PS"))
circos.par(gap.after = c(2, 2, 20))
col_fun <- colorRamp2(c(2, 4), c("#ed8a3f",  "#9b8281"))
circos.heatmap(clinical[,c(1,1)], split = split,
               col = col_fun,track.height = 0.14)

col_fun <- colorRamp2(c(0, 1), c("#d3e2f2",  "#E74C3C"))
circos.heatmap(clinical[,2:5], split = split,
               col = col_fun,track.height = 0.28)
circos.track(
  track.index = get.current.track.index(), 
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 5) {
      cn = colnames(mat)
      n = length(cn)
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
        1:n - 0.5, cn, cex = 0.5, adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }, 
  bg.border = NA
)
circos.clear()

###  Figure 1B right ###
library(circlize)
library(ggplot2)
library(ggforce)
library(ggpubr)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,row.names = 1)
data.type <- data.frame(ratio = as.numeric(table(clinical$type)/93),
                        type = c("LCNEC","Combined","SCLC"))
data.type$type <- factor(data.type$type,
                         levels = c("LCNEC","Combined","SCLC"))
data.tusite <- data.frame(ratio = as.numeric(table(clinical$Tumor_site)/93),
                          type = c("left","left-down","left-up",
                                   "right","right-down","right-middle","right-up"))
data.tusite$type <- factor(data.tusite$type,
                           levels = c("left","left-down","left-up",
                                      "right","right-down","right-middle","right-up"))
data.sex <- data.frame(ratio = as.numeric(table(clinical$sex)/93),
                       type = c("female","male"))
data.sex$type <- factor(data.sex$type,
                        levels = c("female","male"))
data.age <- data.frame(ratio = as.numeric(table(clinical$age_60)/93),
                       type = c("<60",">=60"))
data.age$type <- factor(data.age$type,
                        levels = c("<60",">=60"))
data.smoke <- data.frame(ratio = as.numeric(table(clinical$smoking)/93),
                         type = c("No","Yes"))
data.smoke$type <- factor(data.smoke$type,
                          levels = c("No","Yes"))
data.t <- data.frame(ratio = as.numeric(table(clinical$T_stage)/93),
                     type = c("1","2","3","4"))
data.t$type <- factor(data.t$type,
                      levels = c("1","2","3","4"))
data.n <- data.frame(ratio = as.numeric(table(clinical$N_metastasis)/93),
                     type = c("No","Yes"))
data.n$type <- factor(data.n$type,
                      levels = c("No","Yes"))
data.thera <- data.frame(ratio = as.numeric(table(clinical$Postoperative_adjuvant_chemotherapy)/93),
                         type = c("No","Yes"))
data.thera$type <- factor(data.thera$type,
                          levels = c("No","Yes"))

a1 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#8a9eb5','#ed8a3f','#9b8281'))+
  geom_arc_bar(data=data.type,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a2 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#8a6eaf','#91c0c0','#3a4a7d',
                               '#ffb6c1','#a7a7f2','#a7f2a7','#ffcf48'))+
  geom_arc_bar(data=data.tusite,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a3 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#d3e2f2','#f79990'))+
  geom_arc_bar(data=data.sex,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a4 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#d3e2f2','#f79990'))+
  geom_arc_bar(data=data.age,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a5 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#d3e2f2','#f79990'))+
  geom_arc_bar(data=data.smoke,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a6 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#dee2d1','#91c0c0','#c29e2e','#647370'))+
  geom_arc_bar(data=data.t,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a7 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#d3e2f2','#f79990'))+
  geom_arc_bar(data=data.n,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
a8 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c('#d3e2f2','#f79990'))+
  geom_arc_bar(data=data.thera,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=ratio,fill=type,alpha = 0.4)
  )
ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,nrow = 2,ncol = 4)



###  Tumor purity Figure 1C  ###  
library(ggplot2)
library(ggpubr)
library(gghalves)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,row.names = 1)
clinical$all <- rep("all",dim(clinical)[1])
ggplot() +
  geom_half_boxplot(data = clinical, 
                    aes(x=all, y=Tumor.purity, 
                        fill = all),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = clinical, 
                  aes(x=all, y=Tumor.purity, 
                      color=all), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#FFB6C1'))+
  scale_color_manual(values = c('#FFB6C1'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='Tumor purity %')
min(clinical$Tumor.purity)
max(clinical$Tumor.purity)
median(clinical$Tumor.purity)
clinical$type <- factor(clinical$type,
                        levels = c("LCNEC","LCNEC+SCLC","SCLC"))
ggplot() +
  geom_half_boxplot(data = clinical, 
                    aes(x=type, y=Tumor.purity, 
                        fill = type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = clinical, 
                  aes(x=type, y=Tumor.purity, 
                      color=type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#526187","#ED8A3F","#9B8281"))+
  scale_color_manual(values = c("#526187","#ED8A3F","#9B8281"))+
  stat_compare_means(data = clinical, 
                     aes(x=type, y=Tumor.purity, 
                         fill=type),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='Tumor purity %')
median(clinical$Tumor.purity[which(clinical$type == "SCLC")])



###  oncoplot  Figure 1D ###
library(maftools)
lcnec = read.maf(maf = "lcnec.maf",clinicalData = "clin.tsv")
ls = read.maf(maf = "lcnec+sclc.maf",clinicalData = "clin.tsv")
sclc = read.maf(maf = "sclc.maf",clinicalData = "clin.tsv")
vc_cols = c("#6069B0","#1C7D6A","#D08927","#1e2b4f","#CE5FB5","#F7CAC9","#FD5D67","#53BDA4")
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
                   'Multi_Hit','Frame_Shift_Ins','In_Frame_Ins','Splice_Site',
                   'In_Frame_Del')
typecolors = c("#526187","#ED8A3F","#9B8281")
names(typecolors) = c("LCNEC", "LCNEC+SCLC", "SCLC")
dfscolors = c("#F79990","#A1DFDB")
names(dfscolors) = c("Yes","No")
oscolors = c("#F79990","#A1DFDB")
names(oscolors) = c("Yes","No")
age_60colors = c("#F79990","#A1DFDB")
names(age_60colors) = c("Yes","No")
sexcolors = c("#F79990","#A1DFDB")
names(sexcolors) = c("male","female")
smokcolors = c("#F79990","#A1DFDB")
names(smokcolors) = c("Yes","No")
tcolors = c("#DEE2D1","#91C0C0","#C29E2E","#647370")
names(tcolors) = c("1","2","3","4")
ncolors = c("#F79990","#A1DFDB")
names(ncolors) = c("Yes","No")
thecolors = c("#F79990","#A1DFDB")
names(thecolors) = c("Yes","No")
the_typecolors = c("#FFAD60","#005792","#E8222D","#A1DFDB","#A39E9E")
names(the_typecolors) = c("NSCLC","SCLC","NSCLC+SCLC","No","Unknown")
fabcolors = list(type = typecolors,
                 DFSstate = dfscolors,
                 Osstate = oscolors,
                 age_60 = age_60colors,
                 sex = sexcolors,
                 smoking = smokcolors,
                 T_stage = tcolors,
                 N_metastasis = ncolors,
                 Postoperative_adjuvant_chemotherapy = thecolors,
                 Postoperative_adjuvant_chemotherapy_type = the_typecolors)
plotmafSummary(maf = lcnec, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE)
plotmafSummary(maf = ls, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE)
plotmafSummary(maf = sclc, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE, 
               titvRaw = FALSE)

###  Figure 1D  ###
library(ggplot2)
mat <- data.frame(type = c("Peptides","Identified proteins",
                           "Comparable proteins"),
                  num = c(81120,9707,6979))
mat$type <- factor(mat$type,
                   levels = c("Peptides","Identified proteins",
                              "Comparable proteins"))
ggplot(mat,aes(x = type, y = num,
               fill = type)) +
  geom_col() +
  scale_fill_manual(values = c('#8B8B8B','#BEBEBE','#D9C0E3')) + 
  scale_y_continuous(expand = c(0,0), trans = scales::pseudo_log_trans(), 
                     breaks = c(0,3,10, 30, 100, 250, 600, 
                                1500,5000,20000,80000),
                     limits = c(0, 82000))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.6),
        legend.position="right") + theme_bw()

###  Figure 1F left  ###
library(ggplot2)
library(ggpubr)
data <- read.csv("6979_Protein_normalize_intensity.csv",
                 stringsAsFactors = F,row.names = 1,check.names = F)
count <- c()
for (i in 1:dim(data)[2]){
  count <- c(count,length(which(data[,i] != 0.00001)))
}
protein_count <- data.frame(count = count,
                            patient = colnames(data))
protein_count <- protein_count[rev(order(protein_count$count)),]
protein_count_pl <- protein_count[which(substr(protein_count$patient,1,2) == "PL"),]
protein_count_ls <- protein_count[which(substr(protein_count$patient,1,2) == "LS"),]
protein_count_ps <- protein_count[which(substr(protein_count$patient,1,2) == "PS"),]
protein_count_pl$rank <- 1:dim(protein_count_pl)[1]
protein_count_ls$rank <- 1:dim(protein_count_ls)[1]
protein_count_ps$rank <- 1:dim(protein_count_ps)[1]
protein_count <- rbind(protein_count_pl,
                       protein_count_ls,
                       protein_count_ps)
protein_count$rank <- 1:dim(protein_count)[1]
protein_count$cancer.type <- substr(protein_count$patient,1,2)
ggplot(protein_count, aes(x=rank, y=count, fill=cancer.type)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c('#ED8A3F','#8A9EB5','#9B8281'))+
  theme_bw() +
  geom_hline(yintercept = median(protein_count$count)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("Protein counts")


###  Figure 1F right  ###
library(ggplot2)
library(ggpubr)
data <- read.csv("6979_Protein_normalize_intensity_log2.csv",
                 stringsAsFactors = F,row.names = 1,check.names = F)
data[data == 0.00001] <- NA
data.pl <- data[,which(substr(colnames(data),1,2) == "PL")]
data.ls <- data[,which(substr(colnames(data),1,2) == "LS")]
data.ps <- data[,which(substr(colnames(data),1,2) == "PS")]
data.pl.value <- as.numeric(as.matrix(data.pl))
data.ls.value <- as.numeric(as.matrix(data.ls))
data.ps.value <- as.numeric(as.matrix(data.ps))
bar.mat <- data.frame(protein = c(data.pl.value,
                                  data.ls.value,
                                  data.ps.value),
                      cancer_type = c(rep("PL",length(data.pl.value)),
                                      rep("LS",length(data.ls.value)),
                                      rep("PS",length(data.ps.value))))
bar.mat$cancer_type <- factor(bar.mat$cancer_type,
                              levels = c("PL","LS","PS"))
ggplot(bar.mat, 
       aes(x = protein,fill = cancer_type)) +
  geom_histogram(bins = 200) +
  scale_fill_manual(values = c('#8A9EB5','#ED8A3F','#9B8281'))















