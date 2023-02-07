\name{find.side.gene}
\alias{find.side.gene}
\title{Identify side genes arounding target genes.}
\description{Identify flanking genes that arounding target genes.}
\usage{find.side.gene(PATH_1, PATH_2, range)}
\arguments{
  \item{PATH_1}{path to target genes' bed file.}
  \item{PATH_2}{path to bed file of genes that you want to explore.}
  \item{range}{The range of side genes arounding target genes you want to identify.}
}
\details{Identify side genes arounding target genes.}
\value{
data.frame with columns of: 
seqnames: chromosome of genes;
start: start point of genes;
end: end point of genes;
ID: gene names;
log2foldchange: the expression level of genes
}
\author{Songchen Yang}
\examples{
\dontrun{
PATH1="path/to/bed/file/"
PATH2="path/to/bed/file/"
result<-find.side.gene(PATH_1 = PATH1, PATH_2 = PATH2, range = 2000)
result
}
}