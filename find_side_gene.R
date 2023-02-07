##' Find side genes
##' 
##' Identify side genes around target genes and their log2foldchanges
##' @param PATH_1 The path to bed file of target genes
##' @param PATH_2 The path to bed file of all genes
##' @param range The range of side genes you want to identify.
##' @return A data frame that including log2foldchanges of target genes and their side genes.
##' @author Songchen Yang
##' @export 
find.side.gene<-function(PATH_1, PATH_2, range){
  if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
  TE<-as_tibble(read.delim(PATH_1, header=FALSE))
  gene<-as_tibble(read.delim(PATH_2, header=FALSE))
  TE.loc<-TE[,1:3]
  gene.loc<-gene[,1:3]
  TE.loc<-apply(TE.loc, 2, as.numeric)
  gene.loc<-apply(gene.loc, 2, as.numeric)
  TE[,1:3]<-TE.loc
  gene[,1:3]<-gene.loc
  rm(gene.loc, TE.loc)
  colnames(TE)<-c("seqname", "start", "end", "X_ID", "X_log2FoldChange")
  colnames(gene)<-c("seqname", "start", "end", "Y_ID", "Y_log2FoldChange")
  TE$s_2000<-TE$start - range
  TE$e_2000<-TE$end + range
  D<-data.frame()
  for (i in 1:nrow(gene)){
    for (j in 1:nrow(TE)){
      if (gene$seqname[i] == TE$seqname[j]){
        if(between(gene$start[i], TE$end[j], TE$e_2000[j])==T){
          M<-gene[i,]
          M[,6:7]<-TE[j,4:5]
          D<-rbind(D,M)
        }
      }
    }
  }
  tmp<-D
  rm(D, M)
  D<-data.frame()
  for (i in 1:nrow(gene)){
    for (j in 1:nrow(TE)){
      if (gene$seqname[i] == TE$seqname[j]){
        if(between(gene$end[i], TE$s_2000[j], TE$start[j])==T){
          M<-gene[i,]
          M[,6:7]<-TE[j,4:5]
          D<-rbind(D,M)
        }
      }
    }
  }
  find_side_gene.result<-rbind(tmp, D)
  rm(D, M, tmp)
  return(find_side_gene.result)}



##' correlation of the side genes and target genes

##' Identify the relationship between the side genes and target genes in expression by using #' linear regression model.
##' @param data The result from find.side.gene()
##' @return a formula: y ~ x, and plot picture
##' @export 
lm_picture<-function(data){
  if (!require("ggpmisc", quietly = TRUE))
    install.packages("ggpmisc")
  p<-ggplot(data, aes(x = X_log2FoldChange, y = Y_log2FoldChange)) + 
    geom_point(size = 1) +
    theme_classic() + 
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0, 'lines'),
      legend.position = "top") + 
    stat_smooth(formula = y ~ x, fill = "blue", method = "lm")+
    stat_fit_deviations(formula = y ~ x, color = "skyblue") + 
    stat_poly_eq(use_label(c("eq", "R2", "P")))
  return(p)
}



##' Comparison between two groups

##' Illustrate the linear relation of expression of two group's side genes and target genes at the same time.
##' @param data1 The result from find.side.gene().
##' @param data12 The result from find.side.gene().
##' @param label1 The lable you want to show your plots .
##' @param label2 The lable you want to show your plots .
##' @return plot picture with smooth line
##' @export 
lm_t_pc<-function(data1, data2, lable1, lable2){
  data1$lable<-lable1
  data2$lable<-lable2
  M<-rbind(l1, l2)
  p<-ggplot(data = M, aes(x = X_log2FoldChange, y = Y_log2FoldChange, fill=lable)) + 
    geom_point(aes(color=lable),size = 1) +
    geom_smooth(aes(color=lable),method = "lm", linewidth=1, 
                formula = y ~ x, se = T,
    ) +
    theme_classic() + 
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = unit(0, 'lines'),
      legend.position = "top"
    ) 
  return(p)
}

