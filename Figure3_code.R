
###  Figure 3A  ###
library(maftools)
lcnec = read.maf(maf = "lcnec.maf",clinicalData = "clin.tsv")
ls = read.maf(maf = "lcnec+sclc.maf",clinicalData = "clin.tsv")
sclc = read.maf(maf = "sclc.maf",clinicalData = "clin.tsv")
vc_cols = c("#6069B0","#ABD9FF","#D08927","#1e2b4f",
            "#CE5FB5","#F7CAC9","#FD5D67","#53BDA4")
genes <- read.csv("drivergene.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
deg <- read.csv("deg.csv",
                stringsAsFactors = F,check.names = F)
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
                   'Multi_Hit','Frame_Shift_Ins','In_Frame_Ins','Splice_Site',
                   'In_Frame_Del')
typecolors = c("#526187","#ED8A3F","#9B8281")
names(typecolors) = c("LCNEC", "LCNEC+SCLC", "SCLC")
gene <- intersect(genes$x,
                  ls@gene.summary$Hugo_Symbol)
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
oncoplot(maf =lcnec, fontSize = 0.45 ,top = 50,
         showTumorSampleBarcodes = T,
         SampleNamefontSize=1.2,
         genes = c(gene,"ADGRB1","MSH4"),
         titleFontSize=1.2,
         legendFontSize=1.2,
         colors = vc_cols,
         clinicalFeatures = c('type'),
         sortByAnnotation = TRUE,
         annotationColor = fabcolors,
         writeMatrix=F)# maf
oncoplot(maf =ls, fontSize = 0.45 ,top = 50,
         showTumorSampleBarcodes = T,
         SampleNamefontSize=0.7,
         genes = c(gene,"ADGRB1","MSH4"),
         titleFontSize=1.2,
         legendFontSize=1.2,
         colors = vc_cols,
         clinicalFeatures = c('type'),
         sortByAnnotation = TRUE,
         annotationColor = fabcolors,
         writeMatrix=F)# maf
oncoplot(maf =sclc, fontSize = 0.45 ,top = 50,
         showTumorSampleBarcodes = T,
         SampleNamefontSize=1.7,
         genes = c(gene,"ADGRB1","MSH4"),
         titleFontSize=1.2,
         legendFontSize=1.2,
         colors = vc_cols,
         clinicalFeatures = c('type'),
         sortByAnnotation = TRUE,
         annotationColor = fabcolors,
         writeMatrix=F)# maf



###  Figure 3B ###
library(ggplot2)
library(ggpubr)
gene <- read.csv("driver_gene.csv",
                 stringsAsFactors = F,row.names = 1)
for (i in 1:dim(gene)[2]){
  weizhi <- which(gene[,i] != "aWt")
  gene[weizhi,i] <- "Mutation"
}
prog.gene.pl <- c('FAM135B','RANBP2','RNF213','WNK2','CREBBP')
prog.gene.ls <- c('HERC2','NIPBL','RB1')
prog.gene.ps <- c('ABL2','ERBB4')
all <- unique(c(prog.gene.pl,prog.gene.ls,prog.gene.ps))
gene.mutant <- gene[,all]
protein.raw <- read.csv("6979_Protein_normalize_intensity.csv",
                        stringsAsFactors = F,row.names = 1,check.names = F)
protein.log <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,row.names = 1,check.names = F)
pro <- c("TP53","BCLAF1","LRP1B","RB1","KMT2D","FAT4","FAT1","NOTCH1","FAT3","RNF213",
         "FAT2","HERC2","CDH10","NBEA","FAM135B","PTPRD","GRIN2A","KMT2C","BIRC6","CR1",
         "CREBBP","POLQ","FBN2","NF1","SMARCA4","ZFHX3","ZNF208","ATM","SETD2","BRCA2",
         "EPHA3","ERBB4","ZNF521","RANBP2","NIPBL","ABL2","ABCB1","MAP2","ZNF814","ARID1A",
         "WNK2","MSH4")
prot_inter <- intersect(pro,protein.raw$`Gene Name`)
protein.raw <- protein.raw[match(prot_inter,protein.raw$`Gene Name`),]
rownames(protein.raw) <- protein.raw$`Gene Name`
protein.raw <- protein.raw[,-c(1,2)]  ##  raw used to FC
protein.log <- protein.log[match(prot_inter,protein.log$`Gene name`),]
rownames(protein.log) <- protein.log$`Gene name`
protein.log <- protein.log[,-c(1,2)]  ##  log2 used to P value
protein.raw.pl <- protein.raw[,which(substr(colnames(protein.raw),1,2) == "PL")]
protein.raw.ls <- protein.raw[,which(substr(colnames(protein.raw),1,2) == "LS")]
protein.raw.ps <- protein.raw[,which(substr(colnames(protein.raw),1,2) == "PS")]
protein.log.pl <- protein.log[,colnames(protein.raw.pl)]
protein.log.ls <- protein.log[,colnames(protein.raw.ls)]
protein.log.ps <- protein.log[,colnames(protein.raw.ps)]
gene.mutant.pl <- gene.mutant[colnames(protein.raw.pl),prog.gene.pl]
gene.mutant.ls <- gene.mutant[colnames(protein.raw.ls),prog.gene.ls]
gene.mutant.ps <- gene.mutant[colnames(protein.raw.ps),prog.gene.ps]
pl.fc <- matrix(,dim(gene.mutant.pl)[2],20)
ls.fc <- matrix(,dim(gene.mutant.ls)[2],20)
ps.fc <- matrix(,dim(gene.mutant.ps)[2],20)
pl.p <- matrix(,dim(gene.mutant.pl)[2],20)
ls.p <- matrix(,dim(gene.mutant.ls)[2],20)
ps.p <- matrix(,dim(gene.mutant.ps)[2],20)
### LCNEC
for (i in 1:dim(gene.mutant.pl)[2]){
  gene <- colnames(gene.mutant.pl)[i]
  for(j in 1:20){
    if(length(unique(gene.mutant.pl[,gene])) > 1){
      a <- wilcox.test(as.numeric(protein.log.pl[j,])~gene.mutant.pl[,gene])
      fc <- mean(as.numeric(protein.raw.pl[j,which(gene.mutant.pl[,gene] == "Mutation")]))/
        mean(as.numeric(protein.raw.pl[j,which(gene.mutant.pl[,gene] == "aWt")]))
      pl.p[i,j] <- a$p.value
      pl.fc[i,j] <- fc
    }
    else {
      pl.p[i,j] <- NA
      pl.fc[i,j] <- NA
    }
  }
}
rownames(pl.p) <- colnames(gene.mutant.pl)
colnames(pl.p) <- rownames(protein.log.pl)
rownames(pl.fc) <- colnames(gene.mutant.pl)
colnames(pl.fc) <- rownames(protein.log.pl)
### LCNEC+SCLC
for (i in 1:dim(gene.mutant.ls)[2]){
  gene <- colnames(gene.mutant.ls)[i]
  for(j in 1:20){
    if(length(unique(gene.mutant.ls[,gene])) > 1){
      a <- wilcox.test(as.numeric(protein.log.ls[j,])~gene.mutant.ls[,gene])
      fc <- mean(as.numeric(protein.raw.ls[j,which(gene.mutant.ls[,gene] == "Mutation")]))/
        mean(as.numeric(protein.raw.ls[j,which(gene.mutant.ls[,gene] == "aWt")]))
      ls.p[i,j] <- a$p.value
      ls.fc[i,j] <- fc
    }
    else {
      ls.p[i,j] <- NA
      ls.fc[i,j] <- NA
    }
  }
}
rownames(ls.p) <- colnames(gene.mutant.ls)
colnames(ls.p) <- rownames(protein.log.ls)
rownames(ls.fc) <- colnames(gene.mutant.ls)
colnames(ls.fc) <- rownames(protein.log.ls)
### SCLC
for (i in 1:dim(gene.mutant.ps)[2]){
  gene <- colnames(gene.mutant.ps)[i]
  for(j in 1:20){
    if(length(unique(gene.mutant.ps[,gene])) > 1){
      a <- wilcox.test(as.numeric(protein.log.ps[j,])~gene.mutant.ps[,gene])
      fc <- mean(as.numeric(protein.raw.ps[j,which(gene.mutant.ps[,gene] == "Mutation")]))/
        mean(as.numeric(protein.raw.ps[j,which(gene.mutant.ps[,gene] == "aWt")]))
      ps.p[i,j] <- a$p.value
      ps.fc[i,j] <- fc
    }
    else {
      ps.p[i,j] <- NA
      ps.fc[i,j] <- NA
    }
  }
}
rownames(ps.p) <- colnames(gene.mutant.ps)
colnames(ps.p) <- rownames(protein.log.ps)
rownames(ps.fc) <- colnames(gene.mutant.ps)
colnames(ps.fc) <- rownames(protein.log.ps)
pl <- data.frame(fc = as.numeric(pl.fc),
                 p = as.numeric(pl.p),
                 pro.gene = c(rep(colnames(pl.fc)[1],dim(pl.fc)[1]),rep(colnames(pl.fc)[2],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[3],dim(pl.fc)[1]),rep(colnames(pl.fc)[4],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[5],dim(pl.fc)[1]),rep(colnames(pl.fc)[6],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[7],dim(pl.fc)[1]),rep(colnames(pl.fc)[8],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[9],dim(pl.fc)[1]),rep(colnames(pl.fc)[10],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[11],dim(pl.fc)[1]),rep(colnames(pl.fc)[12],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[13],dim(pl.fc)[1]),rep(colnames(pl.fc)[14],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[15],dim(pl.fc)[1]),rep(colnames(pl.fc)[16],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[17],dim(pl.fc)[1]),rep(colnames(pl.fc)[18],dim(pl.fc)[1]),
                              rep(colnames(pl.fc)[19],dim(pl.fc)[1]),rep(colnames(pl.fc)[20],dim(pl.fc)[1])),
                 mut.gene = rep(rownames(pl.fc),20))
ls <- data.frame(fc = as.numeric(ls.fc),
                 p = as.numeric(ls.p),
                 pro.gene = c(rep(colnames(ls.fc)[1],dim(ls.fc)[1]),rep(colnames(ls.fc)[2],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[3],dim(ls.fc)[1]),rep(colnames(ls.fc)[4],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[5],dim(ls.fc)[1]),rep(colnames(ls.fc)[6],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[7],dim(ls.fc)[1]),rep(colnames(ls.fc)[8],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[9],dim(ls.fc)[1]),rep(colnames(ls.fc)[10],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[11],dim(ls.fc)[1]),rep(colnames(ls.fc)[12],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[13],dim(ls.fc)[1]),rep(colnames(ls.fc)[14],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[15],dim(ls.fc)[1]),rep(colnames(ls.fc)[16],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[17],dim(ls.fc)[1]),rep(colnames(ls.fc)[18],dim(ls.fc)[1]),
                              rep(colnames(ls.fc)[19],dim(ls.fc)[1]),rep(colnames(ls.fc)[20],dim(ls.fc)[1])),
                 mut.gene = rep(rownames(ls.fc),20))
ps <- data.frame(fc = as.numeric(ps.fc),
                 p = as.numeric(ps.p),
                 pro.gene = c(rep(colnames(ps.fc)[1],dim(ps.fc)[1]),rep(colnames(ps.fc)[2],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[3],dim(ps.fc)[1]),rep(colnames(ps.fc)[4],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[5],dim(ps.fc)[1]),rep(colnames(ps.fc)[6],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[7],dim(ps.fc)[1]),rep(colnames(ps.fc)[8],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[9],dim(ps.fc)[1]),rep(colnames(ps.fc)[10],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[11],dim(ps.fc)[1]),rep(colnames(ps.fc)[12],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[13],dim(ps.fc)[1]),rep(colnames(ps.fc)[14],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[15],dim(ps.fc)[1]),rep(colnames(ps.fc)[16],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[17],dim(ps.fc)[1]),rep(colnames(ps.fc)[18],dim(ps.fc)[1]),
                              rep(colnames(ps.fc)[19],dim(ps.fc)[1]),rep(colnames(ps.fc)[20],dim(ps.fc)[1])),
                 mut.gene = rep(rownames(ps.fc),20))
pl$fc[which(is.na(pl$fc) == T)] <- 1
pl$p[which(is.na(pl$p) == T)] <- 1
ls$fc[which(is.na(ls$fc) == T)] <- 1
ls$p[which(is.na(ls$p) == T)] <- 1
ps$fc[which(is.na(ps$fc) == T)] <- 1
ps$p[which(is.na(ps$p) == T)] <- 1
all <- rbind(pl,ls,ps)
all$shape <- rep("tans",dim(all)[1])
all$shape[which(all$pro.gene == all$mut.gene)] <- "cis"

ggplot(data = all, aes(x = mut.gene, y = pro.gene,shape = shape)) +  # 
  geom_point(aes(size = (-log10(p)),color = fc), alpha = 1)+
  scale_shape_manual(values = c(19,15))+
  scale_color_gradient2(low = "#2F3A59",
                        mid = "white",
                        high = "#FF0000",
                        midpoint = 1)


###  Figure 3C  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
gene <- read.csv("driver_gene.csv",
                 stringsAsFactors = F,row.names = 1)
for (i in 1:dim(gene)[2]){
  weizhi <- which(gene[,i] != "aWt")
  gene[weizhi,i] <- "Mutation"
}
pl.gene <- c('CREBBP','RANBP2','RNF213','WNK2')
ls.gene <- c('HERC2','NIPBL','RB1',"NOTCH1")
ps.gene <- c('ABL2')
prot_inter <- c(pl.gene,ls.gene,ps.gene)
protein.log <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,row.names = 1,check.names = F)
protein.log <- protein.log[match(prot_inter,protein.log$`Gene name`),]
rownames(protein.log) <- protein.log$`Gene name`
protein.log <- protein.log[,-c(1,2)]  ##  log2 used to P value
protein.log <- as.data.frame(t(protein.log))
protein.log <- protein.log[,c(pl.gene,ls.gene,ps.gene)]
gene <- gene[,c(pl.gene,ls.gene,ps.gene)]
protein.log.pl <- protein.log[which(substr(rownames(protein.log),1,2) == "PL"),
                              pl.gene]
protein.log.ls <- protein.log[which(substr(rownames(protein.log),1,2) == "LS"),
                              ls.gene]
protein.log.ps <- protein.log[which(substr(rownames(protein.log),1,2) == "PS"),]
gene.pl <- gene[rownames(protein.log.pl),pl.gene]
gene.ls <- gene[rownames(protein.log.ls),ls.gene]
gene.ps <- gene[rownames(protein.log.ps),c(ls.gene,ps.gene)]
protein.log.pl.gene <- cbind(protein.log.pl,gene.pl)
protein.log.ls.gene <- cbind(protein.log.ls,gene.ls)
protein.log.ps.gene <- cbind(protein.log.ps,gene.ps)
colnames(protein.log.pl.gene)[5:8] <- paste(colnames(protein.log.pl.gene)[5:8],
                                            ".label",sep = "")
colnames(protein.log.ls.gene)[5:8] <- paste(colnames(protein.log.ls.gene)[5:8],
                                            ".label",sep = "")
colnames(protein.log.ps.gene)[10:14] <- paste(colnames(protein.log.ps.gene)[10:14],
                                              ".label",sep = "")
pl1 <- ggplot() +
  geom_half_boxplot(data = protein.log.pl.gene, 
                    aes(x=CREBBP.label, y=CREBBP, fill = CREBBP.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.pl.gene, 
                  aes(x=CREBBP.label, y=CREBBP, color=CREBBP.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.pl.gene, 
                     aes(x=CREBBP.label, y=CREBBP, fill=CREBBP.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="CREBBP")
pl2 <- ggplot() +
  geom_half_boxplot(data = protein.log.pl.gene, 
                    aes(x=RANBP2.label, y=RANBP2, fill = RANBP2.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.pl.gene, 
                  aes(x=RANBP2.label, y=RANBP2, color=RANBP2.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.pl.gene, 
                     aes(x=RANBP2.label, y=RANBP2, fill=RANBP2.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="RANBP2")
pl3 <- ggplot() +
  geom_half_boxplot(data = protein.log.pl.gene, 
                    aes(x=RNF213.label, y=RNF213, fill = RNF213.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.pl.gene, 
                  aes(x=RNF213.label, y=RNF213, color=RNF213.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.pl.gene, 
                     aes(x=RNF213.label, y=RNF213, fill=RNF213.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="RNF213")
ls1 <- ggplot() +
  geom_half_boxplot(data = protein.log.ls.gene, 
                    aes(x=HERC2.label, y=HERC2, fill = HERC2.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.ls.gene, 
                  aes(x=HERC2.label, y=HERC2, color=HERC2.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.ls.gene, 
                     aes(x=HERC2.label, y=HERC2, fill=HERC2.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="HERC2")
ls2 <- ggplot() +
  geom_half_boxplot(data = protein.log.ls.gene, 
                    aes(x=NIPBL.label, y=NIPBL, fill = NIPBL.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.ls.gene, 
                  aes(x=NIPBL.label, y=NIPBL, color=NIPBL.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.ls.gene, 
                     aes(x=NIPBL.label, y=NIPBL, fill=NIPBL.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="NIPBL")
ls3 <- ggplot() +
  geom_half_boxplot(data = protein.log.ls.gene, 
                    aes(x=RB1.label, y=RB1, fill = RB1.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.ls.gene, 
                  aes(x=RB1.label, y=RB1, color=RB1.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.ls.gene, 
                     aes(x=RB1.label, y=RB1, fill=RB1.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="RB1")
ls4 <- ggplot() +
  geom_half_boxplot(data = protein.log.ls.gene, 
                    aes(x=RB1.label, y=NOTCH1, fill = RB1.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.ls.gene, 
                  aes(x=RB1.label, y=NOTCH1, color=RB1.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.ls.gene, 
                     aes(x=RB1.label, y=NOTCH1, fill=RB1.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="NOTCH1")
protein.log.ps.gene
ps1 <- ggplot() +
  geom_half_boxplot(data = protein.log.ps.gene, 
                    aes(x=ABL2.label, y=ABL2, fill = ABL2.label),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = protein.log.ps.gene, 
                  aes(x=ABL2.label, y=ABL2, color=ABL2.label), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 2.8)+
  scale_fill_manual(values = c('#CFD3CE','#FF6361'))+
  scale_color_manual(values = c('#CFD3CE','#FF6361'))+
  stat_compare_means(data = protein.log.ps.gene, 
                     aes(x=ABL2.label, y=ABL2, fill=ABL2.label),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="ABL2")
ggarrange(pl1,pl2,pl3,ls1,ls2,ls3,ls4,ps1,nrow = 1,ncol = 8)


###  Figure 3D  ##
library(ggplot2)
library(ggpubr)
library(gghalves)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
gene <- read.csv("driver_gene.csv",
                 stringsAsFactors = F,row.names = 1)
for (i in 1:dim(gene)[2]){
  weizhi <- which(gene[,i] != "aWt")
  gene[weizhi,i] <- "Mutation"
}
data.origin <- read.csv("6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
gene1 <- intersect(colnames(gene),data.origin$`Gene name`) 
data.gene <- data.origin[match(gene1,data.origin$`Gene name`),]
rownames(data.gene) <- data.gene$`Gene name`
data.gene <- data.gene[,-c(1,2)]  
data.gene <- as.data.frame(t(data.gene))
gene <- gene[rownames(data.gene),]
data.gene <- data.gene[rownames(gene),]
colnames(gene) <- paste(colnames(gene),".label",sep = "")
gene.pl <- gene[which(substr(rownames(gene),1,2) == "PL"),]
gene.ls <- gene[which(substr(rownames(gene),1,2) == "LS"),]
gene.ps <- gene[which(substr(rownames(gene),1,2) == "PS"),]
clinical.pl <- clinical[rownames(gene.pl),]
clinical.ls <- clinical[rownames(gene.ls),]
clinical.ps <- clinical[rownames(gene.ps),]
clinical.pl <- cbind(clinical.pl,gene.pl)
clinical.ls <- cbind(clinical.ls,gene.ls)
clinical.ps <- cbind(clinical.ps,gene.ps)
clinical.pl <- clinical.pl[which(clinical.pl$Postoperative_adjuvant_chemotherapy == "Yes"),]
clinical.ls <- clinical.ls[which(clinical.ls$Postoperative_adjuvant_chemotherapy == "Yes"),]
clinical.ps <- clinical.ps[which(clinical.ps$Postoperative_adjuvant_chemotherapy == "Yes"),]
### C-SCLC
smar.t <- table(clinical.ls$RB1.label,clinical.ls$T_stage)
smar.t
fisher.test(smar.t)
smar.n <- table(clinical.ls$RB1.label,clinical.ls$N_metastasis)
smar.n
fisher.test(smar.n)
smar.smok <- table(clinical.ls$RB1.label,clinical.ls$smoking)
smar.smok <- cbind(smar.smok,c(0,0))
fisher.test(smar.smok)


###  Figure 3E left  ###
library(survival)
library(survminer)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
gene <- read.csv("driver_gene.csv",
                 stringsAsFactors = F,row.names = 1)
for (i in 1:dim(gene)[2]){
  weizhi <- which(gene[,i] != "aWt")
  gene[weizhi,i] <- "Mutation"
}
colnames(gene) <- paste(colnames(gene),".label",sep = "")
gene.pl <- gene[which(substr(rownames(gene),1,2) == "PL"),]
gene.ls <- gene[which(substr(rownames(gene),1,2) == "LS"),]
gene.ps <- gene[which(substr(rownames(gene),1,2) == "PS"),]
clinical.pl <- clinical[rownames(gene.pl),]
clinical.ls <- clinical[rownames(gene.ls),]
clinical.ps <- clinical[rownames(gene.ps),]
clinical.pl <- cbind(clinical.pl,gene.pl)
clinical.ls <- cbind(clinical.ls,gene.ls)
clinical.ps <- cbind(clinical.ps,gene.ps)
clinical.pl <- clinical.pl[which(clinical.pl$Postoperative_adjuvant_chemotherapy == "Yes"),]
clinical.ls <- clinical.ls[which(clinical.ls$Postoperative_adjuvant_chemotherapy == "Yes"),]
clinical.ps <- clinical.ps[which(clinical.ps$Postoperative_adjuvant_chemotherapy == "Yes"),]
clinical.pl$SMARCA4.label.new <- paste(clinical.pl$SMARCA4.label,
                                       clinical.pl$Postoperative_adjuvant_chemotherapy)
clinical.pl$RANBP2.label.new <- paste(clinical.pl$RANBP2.label,
                                      clinical.pl$Postoperative_adjuvant_chemotherapy)
clinical.ls$RB1.label.new <- paste(clinical.ls$RB1.label,
                                   clinical.ls$Postoperative_adjuvant_chemotherapy)
clinical.ls$WNK2.label.new <- paste(clinical.ls$WNK2.label,
                                    clinical.ls$Postoperative_adjuvant_chemotherapy)
ls.dfs <- survfit(Surv(DFSm,DFSstate)~RB1.label.new,
                  data = clinical.ls)
summary(ls.dfs)
ggsurvplot(ls.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#CFD3CE","#FF6361"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("WT yes","Mutant yes"), 
           xlab = "Time (months)",ylab = "Disease-free Survival")
table(paste(clinical.ls$DFSstate,clinical.ls$RB1.label))


###  Figure 3E right   ###
library(maxstat)
library(survival)
library(survminer)
clinical <- read.csv("survival data.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
gene1 <- intersect(c("RB1","NOTCH1"),data.origin$`Gene name`) 
data.gene <- data.origin[match(gene1,data.origin$`Gene name`),]
rownames(data.gene) <- data.gene$`Gene name`
data.gene <- data.gene[,-c(1,2)]  
data.gene <- as.data.frame(t(data.gene))
data.gene.ls <- data.gene[which(substr(rownames(data.gene),1,2) == "LS"),]
clinical <- clinical[rownames(data.gene.ls),]
clinical$RB1 <- data.gene.ls$RB1
clinical$NOTCH1 <- data.gene.ls$NOTCH1
clinical$RB1.label <- clinical$RB1
clinical$NOTCH1.label <- clinical$NOTCH1
rb1 <- maxstat.test(Surv(DFSm,DFSstate) ~ RB1, 
                    data = clinical, smethod = "LogRank")
clinical$RB1.label[clinical$RB1 <= rb1$estimate] <- "Low RB1"
clinical$RB1.label[clinical$RB1 > rb1$estimate] <- "High RB1"
RB1.dfs <- survfit(Surv(DFSm,DFSstate)~RB1.label,
                   data = clinical)
summary(RB1.dfs)
ggsurvplot(RB1.dfs,
           pval = TRUE,risk.table = TRUE, risk.table.col = "strata", 
           palette = c("#FF6361","#CFD3CE"),legend = c(2,0.5), 
           legend.title = "", conf.int = F,
           legend.labs = c("High RB1","Low RB1"), 
           xlab = "Time (months)",ylab = "Disease-free Survival")


###  Figure 3F   ###
library(ggplot2)
library(ggpubr)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
data.origin <- read.csv("6979_Protein_normalize_intensity.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
data.ls <- data.origin[,which(substr(colnames(data.origin),1,2) == "LS")]
data.ls <- as.matrix(data.ls)
rownames(data.ls) <- data.origin$`Gene Name`
gene <- read.csv("driver_gene.csv",
                 stringsAsFactors = F,row.names = 1)
for (i in 1:dim(gene)[2]){
  weizhi <- which(gene[,i] != "aWt")
  gene[weizhi,i] <- "Mutation"
}
rb1.ls <- gene[which(substr(rownames(gene),1,2) == 'LS'),c("RB1","RB1")]
rb1.ls.mut <- rownames(rb1.ls)[which(rb1.ls$RB1 == "Mutation")]
rb1.ls.wt <- rownames(rb1.ls)[which(rb1.ls$RB1 == "aWt")]
rb1.ls.mut <- intersect(colnames(data.ls),rb1.ls.mut)
rb1.ls.wt <- intersect(colnames(data.ls),rb1.ls.wt)
fc <- rowMeans(data.ls[,rb1.ls.mut])/rowMeans(data.ls[,rb1.ls.wt])
fc.name <-fc[rev(order(fc))]
gmt <- read.gmt("278_CLUE_result_RB1.gmt")

hh1.res <- GSEA(fc.name, TERM2GENE=gmt,
                minGSSize = 1,
                maxGSSize = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "fdr",
                eps = 1e-10)
hh1.res <- as.data.frame(hh1.res)
hh1.res.6 <- hh1.res[1:6,]
hh1.res.6 <- hh1.res.6[order(hh1.res.6$NES),]
hh1.res.6$Description <- factor(hh1.res.6$Description,
                                levels = c(hh1.res.6$Description))
hh1.res <- hh1.res[order(hh1.res$NES),]
hh1.res$Description <- factor(hh1.res$Description,
                              levels = c(hh1.res$Description))
hh1.res$type <- factor(c("low","low","low",
                         rep("none",length(hh1.res$Description)-6),
                         "high","high","high"),
                       levels = c("high","none","low"))
ggplot(hh1.res, aes(x=Description, y=NES,fill = type)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  scale_fill_manual(values = c('#FF7472','#F9F3D7','#D1D4D0'))+
  theme_bw() +
  geom_hline(yintercept = 0)+
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
  ylab("NES")




























