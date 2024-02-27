
####   AUC   ####
library(pROC)
library(ggplot2)
library(ggpubr)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
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
lcnec_up.gene <- data.origin[lcnec_up,]
ls_up.gene <- data.origin[ls_up,]
sclc_up.gene <- data.origin[sclc_up,]
rownames(lcnec_up.gene) <- lcnec_up.gene$`Gene name`
rownames(ls_up.gene) <- ls_up.gene$`Gene name`
rownames(sclc_up.gene) <- sclc_up.gene$`Gene name`
lcnec_up.gene <- lcnec_up.gene[,-c(1,2)]
ls_up.gene <- ls_up.gene[,-c(1,2)]
sclc_up.gene <- sclc_up.gene[,-c(1,2)]
lcnec_up.gene <- as.data.frame(t(lcnec_up.gene))
ls_up.gene <- as.data.frame(t(ls_up.gene))
sclc_up.gene <- as.data.frame(t(sclc_up.gene))
lcnec_up.gene$cancer_type <- substr(rownames(lcnec_up.gene),1,2)
ls_up.gene$cancer_type <- substr(rownames(ls_up.gene),1,2)
sclc_up.gene$cancer_type <- substr(rownames(sclc_up.gene),1,2)
lcnec_up.gene.plls <- lcnec_up.gene[which(lcnec_up.gene$cancer_type != "PS"),]
lcnec_up.gene.plps <- lcnec_up.gene[which(lcnec_up.gene$cancer_type != "LS"),]
ls_up.gene.plls <- ls_up.gene[which(ls_up.gene$cancer_type != "PS"),]
ls_up.gene.lsps <- ls_up.gene[which(ls_up.gene$cancer_type != "PL"),]
sclc_up.gene.plps <- sclc_up.gene[which(sclc_up.gene$cancer_type != "LS"),]
sclc_up.gene.lsps <- sclc_up.gene[which(sclc_up.gene$cancer_type != "PL"),]
lcnec_up.gene$cancer_type[which(lcnec_up.gene$cancer_type != "PL")] <- "other"
ls_up.gene$cancer_type[which(ls_up.gene$cancer_type != "LS")] <- "other"
sclc_up.gene$cancer_type[which(sclc_up.gene$cancer_type != "PS")] <- "other"
lcnec_up.gene$cancer_type <- factor(lcnec_up.gene$cancer_type,
                                    levels = c("other","PL"))
ls_up.gene$cancer_type <- factor(ls_up.gene$cancer_type,
                                 levels = c("other","LS"))
sclc_up.gene$cancer_type <- factor(sclc_up.gene$cancer_type,
                                   levels = c("other","PS"))
auc.lcnec <- c()
auc.lcnec.plls <- c()
auc.lcnec.plps <- c()
auc.ls <- c()
auc.ls.plls <- c()
auc.ls.lsps <- c()
auc.sclc <- c()
auc.sclc.plps <- c()
auc.sclc.lsps <- c()
## LCNEC
for (i in 1:61){
  a <- roc(lcnec_up.gene$cancer_type,as.numeric(lcnec_up.gene[,i]))
  auc.lcnec <- c(auc.lcnec,a$auc)
}
names(auc.lcnec) <- colnames(lcnec_up.gene)[1:61]
for (i in 1:61){
  a <- roc(lcnec_up.gene.plls$cancer_type,as.numeric(lcnec_up.gene.plls[,i]))
  auc.lcnec.plls <- c(auc.lcnec.plls,a$auc)
}
names(auc.lcnec.plls) <- colnames(lcnec_up.gene.plls)[1:61]
for (i in 1:61){
  a <- roc(lcnec_up.gene.plps$cancer_type,as.numeric(lcnec_up.gene.plps[,i]))
  auc.lcnec.plps <- c(auc.lcnec.plps,a$auc)
}
names(auc.lcnec.plps) <- colnames(lcnec_up.gene.plps)[1:61]
## lcnec+sclc
for (i in 1:9){
  a <- roc(ls_up.gene$cancer_type,as.numeric(ls_up.gene[,i]))
  auc.ls <- c(auc.ls,a$auc)
}
names(auc.ls) <- colnames(ls_up.gene)[1:9]
for (i in 1:9){
  a <- roc(ls_up.gene.plls$cancer_type,as.numeric(ls_up.gene.plls[,i]))
  auc.ls.plls <- c(auc.ls.plls,a$auc)
}
names(auc.ls.plls) <- colnames(ls_up.gene.plls)[1:9]
for (i in 1:9){
  a <- roc(ls_up.gene.lsps$cancer_type,as.numeric(ls_up.gene.lsps[,i]))
  auc.ls.lsps <- c(auc.ls.lsps,a$auc)
}
names(auc.ls.lsps) <- colnames(ls_up.gene.lsps)[1:9]
### sclc
for (i in 1:126){
  a <- roc(sclc_up.gene$cancer_type,as.numeric(sclc_up.gene[,i]))
  auc.sclc <- c(auc.sclc,a$auc)
}
names(auc.sclc) <- colnames(sclc_up.gene)[1:126]
for (i in 1:126){
  a <- roc(sclc_up.gene.plps$cancer_type,as.numeric(sclc_up.gene.plps[,i]))
  auc.sclc.plps <- c(auc.sclc.plps,a$auc)
}
names(auc.sclc.plps) <- colnames(sclc_up.gene.plps)[1:126]
for (i in 1:126){
  a <- roc(sclc_up.gene.lsps$cancer_type,as.numeric(sclc_up.gene.lsps[,i]))
  auc.sclc.lsps <- c(auc.sclc.lsps,a$auc)
}
names(auc.sclc.lsps) <- colnames(sclc_up.gene.lsps)[1:126]
### Merge
auc.lcnec <- auc.lcnec[rev(order(auc.lcnec))]
auc.lcnec.plls <- auc.lcnec.plls[rev(order(auc.lcnec.plls))]
auc.lcnec.plps <- auc.lcnec.plps[rev(order(auc.lcnec.plps))]
auc.ls <- auc.ls[rev(order(auc.ls))]
auc.ls.plls <- auc.ls.plls[rev(order(auc.ls.plls))]
auc.ls.lsps <- auc.ls.lsps[rev(order(auc.ls.lsps))]
auc.sclc <- auc.sclc[rev(order(auc.sclc))]
auc.sclc.plps <- auc.sclc.plps[rev(order(auc.sclc.plps))]
auc.sclc.lsps <- auc.sclc.lsps[rev(order(auc.sclc.lsps))]
auc.lcnec <- data.frame(auc = as.numeric(auc.lcnec),
                        gene = names(auc.lcnec),
                        rank = 1:length(auc.lcnec))
auc.lcnec.plls <- data.frame(auc = as.numeric(auc.lcnec.plls),
                             gene = names(auc.lcnec.plls),
                             rank = 1:length(auc.lcnec.plls))
auc.lcnec.plps <- data.frame(auc = as.numeric(auc.lcnec.plps),
                             gene = names(auc.lcnec.plps),
                             rank = 1:length(auc.lcnec.plps))
auc.ls <- data.frame(auc = as.numeric(auc.ls),
                     gene = names(auc.ls),
                     rank = 1:length(auc.ls))
auc.ls.plls <- data.frame(auc = as.numeric(auc.ls.plls),
                          gene = names(auc.ls.plls),
                          rank = 1:length(auc.ls.plls))
auc.ls.lsps <- data.frame(auc = as.numeric(auc.ls.lsps),
                          gene = names(auc.ls.lsps),
                          rank = 1:length(auc.ls.lsps))
auc.sclc <- data.frame(auc = as.numeric(auc.sclc),
                       gene = names(auc.sclc),
                       rank = 1:length(auc.sclc))
auc.sclc.plps <- data.frame(auc = as.numeric(auc.sclc.plps),
                            gene = names(auc.sclc.plps),
                            rank = 1:length(auc.sclc.plps))
auc.sclc.lsps <- data.frame(auc = as.numeric(auc.sclc.lsps),
                            gene = names(auc.sclc.lsps),
                            rank = 1:length(auc.sclc.lsps))
auc.lcnec <- auc.lcnec[which(auc.lcnec$auc>0.8),]
auc.lcnec.plls <- auc.lcnec.plls[which(auc.lcnec.plls$auc>0.8),]
auc.lcnec.plps <- auc.lcnec.plps[which(auc.lcnec.plps$auc>0.8),]
auc.ls <- auc.ls[1:3,]
auc.ls.plls <- auc.ls.plls[1:3,]
auc.ls.lsps <- auc.ls.lsps[1:3,]
auc.sclc <- auc.sclc[which(auc.sclc$auc>0.8),]
auc.sclc.plps <- auc.sclc.plps[which(auc.sclc.plps$auc>0.8),]
auc.sclc.lsps <- auc.sclc.lsps[which(auc.sclc.lsps$auc>0.8),]
auc.sclc <- auc.sclc[1:10,]
auc.sclc.plps <- auc.sclc.plps[1:10,]
auc.sclc.lsps <- auc.sclc.lsps[1:10,]
###  LCNEC up gene predicted result
auc.lcnec$gene <- factor(auc.lcnec$gene,
                         levels = auc.lcnec$gene)
auc.lcnec.plls$gene <- factor(auc.lcnec.plls$gene,
                              levels = auc.lcnec.plls$gene)
auc.lcnec.plps$gene <- factor(auc.lcnec.plps$gene,
                              levels = auc.lcnec.plps$gene)
pl1 <- ggplot(auc.lcnec, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("LCNEC vs. Others") + coord_flip()
pl2 <- ggplot(auc.lcnec.plls, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("LCNEC vs. Combined") + coord_flip()
pl3 <- ggplot(auc.lcnec.plps, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("LCNEC vs. SCLC") + coord_flip()
###  LS up gene predicted result
auc.ls$gene <- factor(auc.ls$gene,
                      levels = auc.ls$gene)
auc.ls.plls$gene <- factor(auc.ls.plls$gene,
                           levels = auc.ls.plls$gene)
auc.ls.lsps$gene <- factor(auc.ls.lsps$gene,
                           levels = auc.ls.lsps$gene)
ls1 <- ggplot(auc.ls, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("Combined vs. Others") + coord_flip()
ls2 <- ggplot(auc.ls.plls, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("Combined vs. LCNEC") + coord_flip()
ls3 <- ggplot(auc.ls.lsps, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("Combined vs. SCLC") + coord_flip()
###  SCLC up gene predicted result
auc.sclc$gene <- factor(auc.sclc$gene,
                        levels = auc.sclc$gene)
auc.sclc.plps$gene <- factor(auc.sclc.plps$gene,
                             levels = auc.sclc.plps$gene)
auc.sclc.lsps$gene <- factor(auc.sclc.lsps$gene,
                             levels = auc.sclc.lsps$gene)
ps1 <- ggplot(auc.sclc, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("SCLC vs. Others") + coord_flip()
ps2 <- ggplot(auc.sclc.plps, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("SCLC vs. LCNEC") + coord_flip()
ps3 <- ggplot(auc.sclc.lsps, aes(x=gene, y=auc)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("SCLC vs. Combined") + coord_flip()
ggarrange(pl1,pl2,pl3,nrow = 1,ncol = 3)
ggarrange(ls1,ls2,ls3,nrow = 1,ncol = 3)
ggarrange(ps1,ps2,ps3,nrow = 1,ncol = 3)
roc(lcnec_up.gene$cancer_type,
    rowMeans(lcnec_up.gene[,c('NEU1','SERPINB1','PSAP')]))
roc(sclc_up.gene$cancer_type,
    rowMeans(sclc_up.gene[,c('AKT3','PIN1','OTULINL',"ACYP1")]))

lcnec.gene.11 <- c('NEU1','SERPINB1','PSAP','TPP1','SH3BGRL','CTSA',
                   'GGT3P','PCBD1','HOPX','HFE','LTA4H')
sclc.gene.10 <- c('AKT3','PIN1','OTULINL','ACYP1','TBPL1',
                  'HRAS','KIF5C','VEZF1','TLE5','CBX1')
combined3 <- c('VRK2','OXSM','CGAS','NSUN6')
lcnec_up.gene.11 <- lcnec_up.gene[,c(lcnec.gene.11,"cancer_type")]
sclc_up.gene.10 <- sclc_up.gene[,c(sclc.gene.10,"cancer_type")]
ls_up.gene.3 <- ls_up.gene[,c(combined3,"cancer_type")]
###   LCNEC   ###
auc.lenec.gene11 <- matrix(,,3)
for (a in 1:10){       # 2
  for (b in (a+1):11){
    au <- roc(lcnec_up.gene.11$cancer_type,
              rowMeans(lcnec_up.gene.11[,c(a,b)]))
    auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                              c(au$auc,"two-gene",
                                paste(colnames(lcnec_up.gene.11)[a],
                                      colnames(lcnec_up.gene.11)[b],sep = " -- ")))
  }
}
for (a in 1:9){       # 3
  for (b in (a+1):10){
    for (c in (b+1):11){
      au <- roc(lcnec_up.gene.11$cancer_type,
                rowMeans(lcnec_up.gene.11[,c(a,b,c)]))
      auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                c(au$auc,"three-gene",
                                  paste(colnames(lcnec_up.gene.11)[a],
                                        colnames(lcnec_up.gene.11)[b],
                                        colnames(lcnec_up.gene.11)[c],sep = " -- ")))
    }
  }
}
for (a in 1:8){       # 4
  for (b in (a+1):9){
    for (c in (b+1):10){
      for (d in (c+1):11){
        au <- roc(lcnec_up.gene.11$cancer_type,
                  rowMeans(lcnec_up.gene.11[,c(a,b,c,d)]))
        auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                  c(au$auc,"four-gene",
                                    paste(colnames(lcnec_up.gene.11)[a],
                                          colnames(lcnec_up.gene.11)[b],
                                          colnames(lcnec_up.gene.11)[c],
                                          colnames(lcnec_up.gene.11)[d],
                                          sep = " -- ")))
      }
    }
  }
}
for (a in 1:7){       # 5
  for (b in (a+1):8){
    for (c in (b+1):9){
      for (d in (c+1):10){
        for (e in (d+1):11){
          au <- roc(lcnec_up.gene.11$cancer_type,
                    rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e)]))
          auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                    c(au$auc,"five-gene",
                                      paste(colnames(lcnec_up.gene.11)[a],
                                            colnames(lcnec_up.gene.11)[b],
                                            colnames(lcnec_up.gene.11)[c],
                                            colnames(lcnec_up.gene.11)[d],
                                            colnames(lcnec_up.gene.11)[e],
                                            sep = " -- ")))
        }
      }
    }
  }
}
for (a in 1:6){       # 6
  for (b in (a+1):7){
    for (c in (b+1):8){
      for (d in (c+1):9){
        for (e in (d+1):10){
          for (f in (e+1):11){
            au <- roc(lcnec_up.gene.11$cancer_type,
                      rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e,f)]))
            auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                      c(au$auc,"six-gene",
                                        paste(colnames(lcnec_up.gene.11)[a],
                                              colnames(lcnec_up.gene.11)[b],
                                              colnames(lcnec_up.gene.11)[c],
                                              colnames(lcnec_up.gene.11)[d],
                                              colnames(lcnec_up.gene.11)[e],
                                              colnames(lcnec_up.gene.11)[f],
                                              sep = " -- ")))
          }
        }
      }
    }
  }
}
for (a in 1:5){       # 7
  for (b in (a+1):6){
    for (c in (b+1):7){
      for (d in (c+1):8){
        for (e in (d+1):9){
          for (f in (e+1):10){
            for (g in (f+1):11){            
              au <- roc(lcnec_up.gene.11$cancer_type,
                        rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e,f,g)]))
              auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                        c(au$auc,"seven-gene",
                                          paste(colnames(lcnec_up.gene.11)[a],
                                                colnames(lcnec_up.gene.11)[b],
                                                colnames(lcnec_up.gene.11)[c],
                                                colnames(lcnec_up.gene.11)[d],
                                                colnames(lcnec_up.gene.11)[e],
                                                colnames(lcnec_up.gene.11)[f],
                                                colnames(lcnec_up.gene.11)[g],
                                                sep = " -- ")))
            }
          }
        }
      }
    }
  }
}
for (a in 1:4){       # 8 
  for (b in (a+1):5){
    for (c in (b+1):6){
      for (d in (c+1):7){
        for (e in (d+1):8){
          for (f in (e+1):9){
            for (g in (f+1):10){    
              for (h in (g+1):11){
                au <- roc(lcnec_up.gene.11$cancer_type,
                          rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e,f,g,h)]))
                auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                          c(au$auc,"eight-gene",
                                            paste(colnames(lcnec_up.gene.11)[a],
                                                  colnames(lcnec_up.gene.11)[b],
                                                  colnames(lcnec_up.gene.11)[c],
                                                  colnames(lcnec_up.gene.11)[d],
                                                  colnames(lcnec_up.gene.11)[e],
                                                  colnames(lcnec_up.gene.11)[f],
                                                  colnames(lcnec_up.gene.11)[g],
                                                  colnames(lcnec_up.gene.11)[h],
                                                  sep = " -- ")))
              }
            }
          }
        }
      }
    }
  }
}
for (a in 1:3){       # 9 
  for (b in (a+1):4){
    for (c in (b+1):5){
      for (d in (c+1):6){
        for (e in (d+1):7){
          for (f in (e+1):8){
            for (g in (f+1):9){    
              for (h in (g+1):10){
                for (i in (h+1):11){
                  au <- roc(lcnec_up.gene.11$cancer_type,
                            rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e,f,g,h,i)]))
                  auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                            c(au$auc,"nine-gene",
                                              paste(colnames(lcnec_up.gene.11)[a],
                                                    colnames(lcnec_up.gene.11)[b],
                                                    colnames(lcnec_up.gene.11)[c],
                                                    colnames(lcnec_up.gene.11)[d],
                                                    colnames(lcnec_up.gene.11)[e],
                                                    colnames(lcnec_up.gene.11)[f],
                                                    colnames(lcnec_up.gene.11)[g],
                                                    colnames(lcnec_up.gene.11)[h],
                                                    colnames(lcnec_up.gene.11)[i],
                                                    sep = " -- ")))
                }
              }
            }
          }
        }
      }
    }
  }
}
for (a in 1:2){       # 10 
  for (b in (a+1):3){
    for (c in (b+1):4){
      for (d in (c+1):5){
        for (e in (d+1):6){
          for (f in (e+1):7){
            for (g in (f+1):8){    
              for (h in (g+1):9){
                for (i in (h+1):10){
                  for (j in (i+1):11){
                    au <- roc(lcnec_up.gene.11$cancer_type,
                              rowMeans(lcnec_up.gene.11[,c(a,b,c,d,e,f,g,h,i,j)]))
                    auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                                              c(au$auc,"ten-gene",
                                                paste(colnames(lcnec_up.gene.11)[a],
                                                      colnames(lcnec_up.gene.11)[b],
                                                      colnames(lcnec_up.gene.11)[c],
                                                      colnames(lcnec_up.gene.11)[d],
                                                      colnames(lcnec_up.gene.11)[e],
                                                      colnames(lcnec_up.gene.11)[f],
                                                      colnames(lcnec_up.gene.11)[g],
                                                      colnames(lcnec_up.gene.11)[h],
                                                      colnames(lcnec_up.gene.11)[i],
                                                      colnames(lcnec_up.gene.11)[j],
                                                      sep = " -- ")))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
au <- roc(lcnec_up.gene.11$cancer_type,    # 11
          rowMeans(lcnec_up.gene.11[,1:11]))  
auc.lenec.gene11 <- rbind(auc.lenec.gene11,
                          c(au$auc,"eleven-gene",
                            paste(colnames(lcnec_up.gene.11)[1],
                                  colnames(lcnec_up.gene.11)[2],
                                  colnames(lcnec_up.gene.11)[3],
                                  colnames(lcnec_up.gene.11)[4],
                                  colnames(lcnec_up.gene.11)[5],
                                  colnames(lcnec_up.gene.11)[6],
                                  colnames(lcnec_up.gene.11)[7],
                                  colnames(lcnec_up.gene.11)[8],
                                  colnames(lcnec_up.gene.11)[9],
                                  colnames(lcnec_up.gene.11)[10],
                                  colnames(lcnec_up.gene.11)[11],
                                  sep = " -- ")))
auc.lenec.gene11 <- auc.lenec.gene11[-1,]
auc.lenec.gene11 <- as.data.frame(auc.lenec.gene11)
colnames(auc.lenec.gene11) <- c("AUC","Number of gene","Gene")
auc.lenec.gene11$AUC <- as.numeric(auc.lenec.gene11$AUC)
auc.lenec.gene11[which(auc.lenec.gene11$AUC == max(auc.lenec.gene11$AUC)),]
###  Combined  ###
au1 <- roc(ls_up.gene.3$cancer_type,    # 2 
          rowMeans(ls_up.gene.3[,1:2])) 
au2 <- roc(ls_up.gene.3$cancer_type,    # 2 
          rowMeans(ls_up.gene.3[,c(1,3)])) 
au3 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(1,4)])) 
au4 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(2,3)])) 
au5 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(2,4)])) 
au6 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(3,4)])) 
au7 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(1,2,3)])) 
au8 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(1,2,4)])) 
au9 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(2,3,4)])) 
au10 <- roc(ls_up.gene.3$cancer_type,    # 2 
           rowMeans(ls_up.gene.3[,c(1,2,3,4)])) 
auc.ls.gene.3 <- data.frame(AUC = c(au1$auc,au2$auc,au3$auc,au4$auc,
                                    au5$auc,au6$auc,au7$auc,au8$auc,
                                    au9$auc,au10$auc),
                            Gene = c("VRK2--OXSM","VRK2--CGAS","VRK2--NSUN6",
                                     "OXSM--CGAS","OXSM--NSUN6","CGAS--NSUN6",
                                     "VRK2--OXSM--CGAS","VRK2--OXSM--NSUN6",
                                     "OXSM--CGAS--NSUN6","VRK2--OXSM--CGAS--NSUN6"))
###  SCLC   ###
auc.sclc.gene10 <- matrix(,,3)
for (a in 1:9){       # 2
  for (b in (a+1):10){
    au <- roc(sclc_up.gene.10$cancer_type,
              rowMeans(sclc_up.gene.10[,c(a,b)]))
    auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                             c(au$auc,"two-gene",
                               paste(colnames(sclc_up.gene.10)[a],
                                     colnames(sclc_up.gene.10)[b],sep = " -- ")))
  }
}
for (a in 1:8){       # 3
  for (b in (a+1):9){
    for (c in (b+1):10){
      au <- roc(sclc_up.gene.10$cancer_type,
                rowMeans(sclc_up.gene.10[,c(a,b,c)]))
      auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                               c(au$auc,"three-gene",
                                 paste(colnames(sclc_up.gene.10)[a],
                                       colnames(sclc_up.gene.10)[b],
                                       colnames(sclc_up.gene.10)[c],
                                       sep = " -- ")))
    }
  }
}
for (a in 1:7){       # 4
  for (b in (a+1):8){
    for (c in (b+1):9){
      for (d in (c+1):10){
        au <- roc(sclc_up.gene.10$cancer_type,
                  rowMeans(sclc_up.gene.10[,c(a,b,c,d)]))
        auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                 c(au$auc,"four-gene",
                                   paste(colnames(sclc_up.gene.10)[a],
                                         colnames(sclc_up.gene.10)[b],
                                         colnames(sclc_up.gene.10)[c],
                                         colnames(sclc_up.gene.10)[d],
                                         sep = " -- ")))
      }
    }
  }
}
for (a in 1:6){       # 5
  for (b in (a+1):7){
    for (c in (b+1):8){
      for (d in (c+1):9){
        for (e in (d+1):10){
          au <- roc(sclc_up.gene.10$cancer_type,
                    rowMeans(sclc_up.gene.10[,c(a,b,c,d,e)]))
          auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                   c(au$auc,"five-gene",
                                     paste(colnames(sclc_up.gene.10)[a],
                                           colnames(sclc_up.gene.10)[b],
                                           colnames(sclc_up.gene.10)[c],
                                           colnames(sclc_up.gene.10)[d],
                                           colnames(sclc_up.gene.10)[e],
                                           sep = " -- ")))
        }
      }
    }
  }
}
for (a in 1:5){       # 6
  for (b in (a+1):6){
    for (c in (b+1):7){
      for (d in (c+1):8){
        for (e in (d+1):9){
          for (f in (e+1):10){
            au <- roc(sclc_up.gene.10$cancer_type,
                      rowMeans(sclc_up.gene.10[,c(a,b,c,d,e,f)]))
            auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                     c(au$auc,"six-gene",
                                       paste(colnames(sclc_up.gene.10)[a],
                                             colnames(sclc_up.gene.10)[b],
                                             colnames(sclc_up.gene.10)[c],
                                             colnames(sclc_up.gene.10)[d],
                                             colnames(sclc_up.gene.10)[e],
                                             colnames(sclc_up.gene.10)[f],
                                             sep = " -- ")))
          }
        }
      }
    }
  }
}
for (a in 1:4){       # 7 
  for (b in (a+1):5){
    for (c in (b+1):6){
      for (d in (c+1):7){
        for (e in (d+1):8){
          for (f in (e+1):9){
            for (g in (f+1):10){
              au <- roc(sclc_up.gene.10$cancer_type,
                        rowMeans(sclc_up.gene.10[,c(a,b,c,d,e,f,g)]))
              auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                       c(au$auc,"seven-gene",
                                         paste(colnames(sclc_up.gene.10)[a],
                                               colnames(sclc_up.gene.10)[b],
                                               colnames(sclc_up.gene.10)[c],
                                               colnames(sclc_up.gene.10)[d],
                                               colnames(sclc_up.gene.10)[e],
                                               colnames(sclc_up.gene.10)[f],
                                               colnames(sclc_up.gene.10)[g],
                                               sep = " -- ")))
            }
          }
        }
      }
    }
  }
}
for (a in 1:3){       # 8
  for (b in (a+1):4){
    for (c in (b+1):5){
      for (d in (c+1):6){
        for (e in (d+1):7){
          for (f in (e+1):8){
            for (g in (f+1):9){
              for (h in (g+1):10){
                au <- roc(sclc_up.gene.10$cancer_type,
                          rowMeans(sclc_up.gene.10[,c(a,b,c,d,e,f,g,h)]))
                auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                         c(au$auc,"eight-gene",
                                           paste(colnames(sclc_up.gene.10)[a],
                                                 colnames(sclc_up.gene.10)[b],
                                                 colnames(sclc_up.gene.10)[c],
                                                 colnames(sclc_up.gene.10)[d],
                                                 colnames(sclc_up.gene.10)[e],
                                                 colnames(sclc_up.gene.10)[f],
                                                 colnames(sclc_up.gene.10)[g],
                                                 colnames(sclc_up.gene.10)[h],
                                                 sep = " -- ")))
              }
            }
          }
        }
      }
    }
  }
}
for (a in 1:2){       # 9
  for (b in (a+1):3){
    for (c in (b+1):4){
      for (d in (c+1):5){
        for (e in (d+1):6){
          for (f in (e+1):7){
            for (g in (f+1):8){
              for (h in (g+1):9){
                for (i in (h+1):10){
                  au <- roc(sclc_up.gene.10$cancer_type,
                            rowMeans(sclc_up.gene.10[,c(a,b,c,d,e,f,g,h,i)]))
                  auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                                           c(au$auc,"nine-gene",
                                             paste(colnames(sclc_up.gene.10)[a],
                                                   colnames(sclc_up.gene.10)[b],
                                                   colnames(sclc_up.gene.10)[c],
                                                   colnames(sclc_up.gene.10)[d],
                                                   colnames(sclc_up.gene.10)[e],
                                                   colnames(sclc_up.gene.10)[f],
                                                   colnames(sclc_up.gene.10)[g],
                                                   colnames(sclc_up.gene.10)[h],
                                                   colnames(sclc_up.gene.10)[i],
                                                   sep = " -- ")))
                }
              }
            }
          }
        }
      }
    }
  }
}
au <- roc(sclc_up.gene.10$cancer_type,  # 10
          rowMeans(sclc_up.gene.10[,1:10]))
auc.sclc.gene10 <- rbind(auc.sclc.gene10,
                         c(au$auc,"ten-gene",
                           paste(colnames(sclc_up.gene.10)[1],
                                 colnames(sclc_up.gene.10)[2],
                                 colnames(sclc_up.gene.10)[3],
                                 colnames(sclc_up.gene.10)[4],
                                 colnames(sclc_up.gene.10)[5],
                                 colnames(sclc_up.gene.10)[6],
                                 colnames(sclc_up.gene.10)[7],
                                 colnames(sclc_up.gene.10)[8],
                                 colnames(sclc_up.gene.10)[9],
                                 colnames(sclc_up.gene.10)[10],
                                 sep = " -- ")))
auc.sclc.gene10 <- auc.sclc.gene10[-1,]
auc.sclc.gene10 <- as.data.frame(auc.sclc.gene10)
colnames(auc.sclc.gene10) <- c("AUC","Number of gene","Gene")
auc.sclc.gene10$AUC <- as.numeric(auc.sclc.gene10$AUC)
auc.sclc.gene10[which(auc.sclc.gene10$AUC == max(auc.sclc.gene10$AUC)),]

auc.lenec.gene11.1 <- auc.lenec.gene11[order(auc.lenec.gene11$AUC),]
auc.lenec.gene11.1$cou <- 1:dim(auc.lenec.gene11)[1]
auc.lenec.gene11.1$`Number of gene` <- factor(auc.lenec.gene11.1$`Number of gene`,
                                              levels = c("two-gene","three-gene","four-gene",
                                                         "five-gene","six-gene","seven-gene",
                                                         "eight-gene","nine-gene","ten-gene",
                                                         "eleven-gene"))
a1 <- ggplot(data = auc.lenec.gene11.1, 
             mapping = aes(x = cou, 
                           y = AUC)) + 
  geom_line() + 
  geom_point(aes(color = `Number of gene`)) +
  scale_color_manual(values = c('#FFCF48','#F7CAC9',
                                '#F7A072','#ff6361','#A8DADC',
                                '#87CEEB','#2A9D8F','#A7A9AC',
                                '#7A4D7B','#7F8FA6'))

auc.sclc.gene10.1 <- auc.sclc.gene10[order(auc.sclc.gene10$AUC),]
auc.sclc.gene10.1$cou <- 1:dim(auc.sclc.gene10)[1]
auc.sclc.gene10.1$`Number of gene` <- factor(auc.sclc.gene10.1$`Number of gene`,
                                             levels = c("two-gene","three-gene","four-gene",
                                                        "five-gene","six-gene","seven-gene",
                                                        "eight-gene","nine-gene","ten-gene"))
a2 <- ggplot(data = auc.sclc.gene10.1, 
             mapping = aes(x = cou, 
                           y = AUC)) + 
  geom_line() + 
  geom_point(aes(color = `Number of gene`)) +
  scale_color_manual(values = c('#FFCF48','#F7CAC9',
                                '#F7A072','#ff6361','#A8DADC',
                                '#87CEEB','#2A9D8F','#A7A9AC',
                                '#7A4D7B','#7F8FA6'))
ggarrange(a1,a2,nrow = 2,ncol = 1)
auc.lenec.gene11.1.top5 <- auc.lenec.gene11.1[2032:2036,]
auc.sclc.gene10.1.top5 <- auc.sclc.gene10.1[1009:1013,]

ggplot(auc.lenec.gene11.1.top5, 
       aes(x=Gene,y=AUC)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("AUC")
ggplot(auc.sclc.gene10.1.top5, 
       aes(x=Gene,y=AUC)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black",
                                  size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab("AUC")
ggplot(auc.ls.gene.3, 
       aes(x=Gene,y=AUC)) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black", 
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("AUC")
library("pROC")
# PL 
roc1 <- roc(lcnec_up.gene$cancer_type,
            rowMeans(lcnec_up.gene[,c('NEU1','SERPINB1')]),
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 1, .01), boot.n=100)  
roc2 <- roc(lcnec_up.gene$cancer_type,
            lcnec_up.gene[,c('NEU1')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj2 <- ci.sp(roc2, sensitivities=seq(0, 1, .01), boot.n=100) 
roc3 <- roc(lcnec_up.gene$cancer_type,
            lcnec_up.gene[,c('SERPINB1')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj3 <- ci.sp(roc3, sensitivities=seq(0, 1, .01), boot.n=100)  
plot(sp.obj1, type="shape", col="#8A9EB5",ci=TRUE)  
plot(sp.obj2, type="shape", col="#F2F2F2",ci=TRUE)  
plot(sp.obj3, type="shape", col="#F2F2F2",ci=TRUE)  #
roc.test(roc1,roc2, method=c("delong"))
roc.test(roc1,roc3, method=c("delong"))
roc.test(roc2,roc3, method=c("delong"))
# LS
roc1 <- roc(ls_up.gene$cancer_type,
            rowMeans(ls_up.gene[,c('VRK2','OXSM','NSUN6')]),
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 1, .01), boot.n=100)  
roc2 <- roc(ls_up.gene$cancer_type,
            ls_up.gene[,c('VRK2')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj2 <- ci.sp(roc2, sensitivities=seq(0, 1, .01), boot.n=100)  
roc3 <- roc(ls_up.gene$cancer_type,
            ls_up.gene[,c('OXSM')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj3 <- ci.sp(roc3, sensitivities=seq(0, 1, .01), boot.n=100)  
roc4 <- roc(ls_up.gene$cancer_type,
            ls_up.gene[,c('NSUN6')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj4 <- ci.sp(roc4, sensitivities=seq(0, 1, .01), boot.n=100)  
plot(sp.obj1, type="shape", col="#F7CDAC",ci=TRUE)  
plot(sp.obj2, type="shape", col="#F2F2F2",ci=TRUE)  
plot(sp.obj3, type="shape", col="#F2F2F2",ci=TRUE)  
plot(sp.obj4, type="shape", col="#F2F2F2",ci=TRUE)  
roc.test(roc1,roc2, method=c("delong"))
roc.test(roc1,roc3, method=c("delong"))
roc.test(roc1,roc4, method=c("delong"))

# PS
roc1 <- roc(sclc_up.gene$cancer_type,
            rowMeans(sclc_up.gene[,c('AKT3','OTULINL','TBPL1')]),
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 1, .01), boot.n=100)  
roc2 <- roc(sclc_up.gene$cancer_type,
            sclc_up.gene[,c('AKT3')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj2 <- ci.sp(roc2, sensitivities=seq(0, 1, .01), boot.n=100)  
roc3 <- roc(sclc_up.gene$cancer_type,
            sclc_up.gene[,c('OTULINL')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj3 <- ci.sp(roc3, sensitivities=seq(0, 1, .01), boot.n=100)  
roc4 <- roc(sclc_up.gene$cancer_type,
            sclc_up.gene[,c('TBPL1')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj4 <- ci.sp(roc4, sensitivities=seq(0, 1, .01), boot.n=100)  
plot(sp.obj1, type="shape", col="#8A9EB5",ci=TRUE)  
plot(sp.obj2, type="shape", col="#F2F2F2",ci=TRUE)  
plot(sp.obj3, type="shape", col="#F2F2F2",ci=TRUE)  
plot(sp.obj4, type="shape", col="#F2F2F2",ci=TRUE)  
roc.test(roc1,roc2, method=c("delong"))
roc.test(roc1,roc3, method=c("delong"))
roc.test(roc1,roc4, method=c("delong"))
roc.test(roc2,roc3, method=c("delong"))
roc.test(roc2,roc4, method=c("delong"))
roc.test(roc3,roc4, method=c("delong"))

roc(ls_up.gene$cancer_type,
    ls_up.gene[,c('VRK2')],
    plot=TRUE, print.thres=TRUE, ci=TRUE,
    print.auc=TRUE,legacy.axes = TRUE,col = "red")

###  heatmap  ###
library(pheatmap)
clini <- read.csv("survival data.csv",
                  stringsAsFactors = F,
                  row.names = 1,check.names = F)
data.origin <- read.csv("dreamai_6979_Protein_normalize_intensity_log2.csv",
                        stringsAsFactors = F,
                        row.names = 1,check.names = F)
gene <- c('NEU1','SERPINB1',  # LCNEC 
          'AKT3','OTULINL','TBPL1',    # SCLC
          'VRK2','OXSM',"NSUN6")    # C-SCLC
data.auc.gene <- data.origin[match(gene,data.origin$`Gene name`),]
rownames(data.auc.gene) <- data.auc.gene$`Gene name`
data.auc.gene <- data.auc.gene[,-c(1,2)]
mean.mat <- rbind(colMeans(data.auc.gene[1:2,]),
                  colMeans(data.auc.gene[3:5,]),
                  colMeans(data.auc.gene[6:8,]))
rownames(mean.mat) <- c("NEU1 & SERPINB1",
                        "AKT3 & OTULINL & TBPL1",
                        "VRK2 & OXSM & NSUN6")
linshi <- apply(data.auc.gene, 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data.auc.gene)
linshi[linshi>1] <- 1
linshi[linshi< (-1)] <- (-1)
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
yellow <- "#FFF200";
black <- "#1F1F1F";
a <- "#4A4A4A"
b <- "#FFF200"
clini <- clini[colnames(linshi),]
annotation_col <- data.frame(Cancer_type = as.factor(substr(rownames(clini),1,2)),
                             DFSstate = as.factor(clini$DFSstate),
                             Osstate = as.factor(clini$Osstate),
                             age_60 = as.factor(clini$age_60),
                             sex = as.factor(clini$sex),
                             smoking = as.factor(clini$smoking),
                             Tumorsite = as.factor(clini$Tumor_site),
                             T_stage = as.factor(clini$T_stage),
                             N_metastasis = as.factor(clini$N_metastasis),
                             M_metastasis = as.factor(clini$M_metastasis),
                             Postoperative_adjuvant_chemotherapy = as.factor(clini$Postoperative_adjuvant_chemotherapy),
                             Postoperative_adjuvant_chemotherapy_type = as.factor(clini$Postoperative_adjuvant_chemotherapy_type))
rownames(annotation_col) <- rownames(clini)
ann_colors = list(Cancer_type = c("PL" = "#8A9EB5","LS" = "#ED8A3F","PS" = "#9B8281"),
                  DFSstate = c("0" = "#A1DFDB", "1" = "#F79990"),
                  Osstate = c("0" = "#A1DFDB", "1" = "#F79990"),
                  age_60 = c("<60" = "#A1DFDB",">=60" = "#F79990"),
                  sex = c("female" = "#A1DFDB","male" = "#F79990"),
                  smoking = c("No"="#A1DFDB","Yes"="#F79990"),
                  Tumorsite = c("left" = "#8A6EAF","left-down" = "#3A4A7D","left-up" = "#7F8FA6",
                                "right" = "#FFB6C1","right-down" = "#FFCF48","right-middle" = "#a7f2a7",
                                "right-up" = "#a7a7f2"),
                  T_stage = c("1"="#dee2d1","2"="#91c0c0",
                              "3"="#c29e2e","4"="#647370"),
                  N_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  M_metastasis = c("No"="#A1DFDB","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy = c("No"="#A1DFDB","Yes"="#F79990"),
                  Postoperative_adjuvant_chemotherapy_type = c("No"="#A1DFDB","NSCLC" = "#ffad60",
                                                               "SCLC" = "#005792","SCLC+NSCLC" = "#E8222D",
                                                               "Unknown" = "#a39e9e"))

pheatmap(linshi,fontsize=6,gaps_col = c(16,46),
         color  = colorRampPalette(c(blue,white,red))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)
mean.linshi <- apply(mean.mat, 1, scale)
mean.linshi <- t(mean.linshi)
colnames(mean.linshi) <- colnames(mean.mat)
mean.linshi[mean.linshi>1] <- 1
mean.linshi[mean.linshi< (-1)] <- (-1)
pheatmap(mean.linshi,
         fontsize=6,gaps_col = c(16,46),
         color  = colorRampPalette(c(black,yellow))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = F)



################################################  
####           IHC data validation          ####
################################################  
###   data process  ###
data <- read.csv("IHC-data.csv",stringsAsFactors = F,
                 check.names = F)
unique_name <- unique(data$ID)
mat <- matrix(,,8)
for (i in 1:length(unique_name)){
  weizhi <- which(unique_name[i] == data$ID)
  if (length(weizhi > 1)){
    mat <- rbind(mat,c(unique_name[i],
                       colMeans(data[weizhi,2:7]),
                       data[weizhi[1],8]))
  }
}
mat <- mat[-1,]
#write.csv(mat,"IHC-data.csv",quote = F)


##  boxplot
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("IHC-data_zzc.csv",stringsAsFactors = F,
                 check.names = F,row.names = 1)
data$NEU1_SERPINB1 <- rowMeans(data[,c('NEU1','SERPINB1')])
ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Histological_type, y=NEU1_SERPINB1, fill = Histological_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Histological_type, y=NEU1_SERPINB1, color=Histological_type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#E76F51","#AEB6FF"))+
  scale_color_manual(values = c("#E76F51","#AEB6FF"))+
  stat_compare_means(data = data, 
                     aes(x=Histological_type, y=NEU1_SERPINB1),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='NEU1 & SERPINB1 (average)')


##   AUC 
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(modEvA)
library(caret)
data_select <- read.csv("IHC-data.csv",stringsAsFactors = F,
                        check.names = F,row.names = 1)
##  LCNEC
roc1 <- roc(data_select$LCNEC,
            rowMeans(data_select[,c('NEU1','SERPINB1')]),
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj1 <- ci.sp(roc1, sensitivities=seq(0, 1, .01), boot.n=100)  ## 
plot(sp.obj1, type="shape", col="#8A9EB5",ci=TRUE)  ## 
AUC(obs=as.character(data_select$LCNEC),
    pred=as.numeric(rowMeans(data_select[,c('NEU1','SERPINB1')]))/300,
    curve ="PR",simplif=T, main = "PR curve")
AUC(obs=as.character(data_select$LCNEC),
    pred=as.numeric(data_select[,c('NEU1')])/300,
    curve ="PR",simplif=T, main = "PR curve")
AUC(obs=as.character(data_select$LCNEC),
    pred=as.numeric(data_select[,c('SERPINB1')])/300,
    curve ="PR",simplif=T, main = "PR curve")
roc2 <- roc(data_select$LCNEC,
            data_select[,c('NEU1')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj2 <- ci.sp(roc2, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj2, type="shape", col="#F2F2F2",ci=TRUE) 
roc3 <- roc(data_select$LCNEC,
            data_select[,c('SERPINB1')],
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj3 <- ci.sp(roc3, sensitivities=seq(0, 1, .01), boot.n=100)  
plot(sp.obj3, type="shape", col="#F2F2F2",ci=TRUE)  
roc.test(roc1,roc2, method=c("delong"))
roc.test(roc1,roc3, method=c("delong"))
roc.test(roc2,roc3, method=c("delong"))
##  SCLC
roc22 <- roc(data_select$SCLC,
             data_select[,c('AKT3 H-score')],
             plot=TRUE, print.thres=TRUE, ci=TRUE,
             print.auc=TRUE,legacy.axes = TRUE,col = "red")
sp.obj22 <- ci.sp(roc22, sensitivities=seq(0, 1, .01), boot.n=100)  


##  boxplt  
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("IHC-data_zzc.csv",stringsAsFactors = F,
                 check.names = F,row.names = 1)
data <- data[which(data$Histological_type != "C-SCLC"),]
data$NEU1_SERPINB1 <- rowMeans(data[,c('NEU1','SERPINB1')])
a1 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Histological_type, y=NEU1_SERPINB1, fill = Histological_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Histological_type, y=NEU1_SERPINB1, color=Histological_type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#8A9EB5","#9B8281"))+
  scale_color_manual(values = c("#8A9EB5","#9B8281"))+
  stat_compare_means(data = data, 
                     aes(x=Histological_type, y=NEU1_SERPINB1),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='NEU1 & SERPINB1 (average H-score)')
a2 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Histological_type, y=NEU1, fill = Histological_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Histological_type, y=NEU1, color=Histological_type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#8A9EB5","#9B8281"))+
  scale_color_manual(values = c("#8A9EB5","#9B8281"))+
  stat_compare_means(data = data, 
                     aes(x=Histological_type, y=NEU1),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='NEU1 H-score')
a3 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Histological_type, y=SERPINB1, fill = Histological_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Histological_type, y=SERPINB1, color=Histological_type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#8A9EB5","#9B8281"))+
  scale_color_manual(values = c("#8A9EB5","#9B8281"))+
  stat_compare_means(data = data, 
                     aes(x=Histological_type, y=SERPINB1),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='SERPINB1 H-score')
a4 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=Histological_type, y=`AKT3 H-score`, fill = Histological_type),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=Histological_type, y=`AKT3 H-score`, color=Histological_type), 
                  transformation = position_jitter(width = 0.1,height = 0),
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c("#8A9EB5","#9B8281"))+
  scale_color_manual(values = c("#8A9EB5","#9B8281"))+
  stat_compare_means(data = data, 
                     aes(x=Histological_type, y=`AKT3 H-score`),
                     transformation = position_jitter(width = 0.1,height = 0),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="right")+
  labs(y='AKT3 H-score')
ggarrange(a1,a2,a3,a4,nrow = 2,ncol = 2)









