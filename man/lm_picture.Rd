\name{lm_picture}
\alias{lm_picture}
\title{correlation between expression levels of target genes and their flanking genes}
\description{lm_picture demonstrate correlations between expression level of target genes  and flanking genes by using linear regression model.}
\usage{lm_picture(data)}
\arguments{
  \item{data}{The data.frame generated from find.side.gene()}
}
\details{demonstration of the relationship between the expression levels of target genes and their side genes by using linear regression model.}
\value{
data.frame with columns of:
log2foldchange of X;
log2foldchange of Y
}
\author{Songchen Yang}
\examples{
\dontrun{
result<-find.side.gene(Path1, Path2, 2000)
lm_picture(result)
}
}

