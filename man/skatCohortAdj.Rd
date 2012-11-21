\name{skatCohortAdj}
\alias{skatCohortAdj}
\alias{skatFamCohortAdj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run SKAT on data from a single cohort, conditional on specified SNP effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function works exactly as \code{\link{skatCohort}}, but with the additional argument `adjustments' specifying genes for which conditional analyses are desired, and which SNPs to condition on.
}
\usage{
skatCohortAdj(Z, formula, family = gaussian(), SNPInfo=NULL, adjustments= NULL, 
	snpNames = "Name", aggregateBy = "gene", data=parent.frame())

skatFamCohortAdj(Z, formula, SNPInfo=NULL, adjustments= NULL, snpNames = "Name", 
	aggregateBy = "gene",fullkins, sparse = TRUE, data=parent.frame())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
	A genotype matrix (dosage matrix) - rows correspond to individuals and columns correspond to SNPs. Use 'NA' for missing values. The column names of this matrix should correspond to SNP names in the SNP information file.
	}
  \item{formula}{
	Base formula, of the kind used in glm() - typically of the form y~covariate1 + covariate2. For Cox models, the formula follows that of the coxph() function.
}
  \item{family}{
	for skatCohort: either gaussian(), for continuous data, or binomial() for 0/1 outcomes. Binary outcomes are not currently supported for family data.
	}
  \item{SNPInfo}{
	SNP Info file - must contain fields given in 'snpName' and 'aggregateBy'. 	
	}
  \item{adjustments}{
	A data frame of the same format at SNPInfo, pairing genes to analyze with snp 
  }
  \item{snpNames}{
	The field of SNPInfo where the SNP identifiers are found.  Default is 'Name'
	}
  \item{aggregateBy}{
    The field of SNPInfo on which the skat results were aggregated.  Default is 'gene'. For single snps which are intended only for single variant analyses, it is reccomended that they have a unique identifier in this field.  
  }
  \item{data}{
	data frame in which to find variables in the formula
  }	
  \item{fullkins}{
  	for skatFamCohort: the kinship matrix. See lmekin and the kinship2 package for more details
  }
  \item{sparse}{
  	for skatFamCohort: whether or not to use a sparse Matrix approximation for dense kinship matrices (defaults to TRUE)
  }
}
\details{
	This function has the same syntax as \code{\link{skatCohort}} and \code{\link{skatFamCohort}}, but requires an extra argument `adjustments`. This is a data frame of the same format as SKAT meta, i.e. with a `snpNames` and `aggregateBy` columns. The function works by looping through the genes in the adjustment file, adding the corresponding SNPs to the null model in SKAT.  For instance, if one wants to adjuste `gene1` for SNPs a and b (which need not be in gene 1), and `gene2' for SNPs c, the adjustments would be something like
	\code{adjustments = data.frame(Name = c("a","b","c"), gene = c("gene1","gene1","gene2"))}

See the examples for an illustration.
}
\value{
	an object of class 'skatCohort'. Note that unlike output from the function \code{\link{skatCohort}}, the null models in each element of the list may be different. When meta analyzing these, it may be good to subset the SNPInfo file to the genes of interest.
}
\author{
Arie Voorman, Jennifer Brody
}

\seealso{
\code{\link{skatCohort}}
\code{\link{skatMeta}}
\code{\link{burdenMeta}}
\code{\link{singlesnpMeta}}
\code{\link{coxph}}
}
\examples{
	
###load example data for two cohorts:
### see ?skatExample	
data(skatExample)

#specify adjustment variables
adjustments <- SNPInfo[c(1:3, 20,100), ]
adjustments


####run on each cohort:
cohort1.adj <- skatCohortAdj(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, 
	adjustments=adjustments, data =pheno1)
cohort2.adj <- skatFamCohortAdj(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
	adjustments=adjustments, fullkins=kins, data=pheno2)

SNPInfo.sub <- subset(SNPInfo, (SNPInfo$gene \%in\% adjustments$gene) & 
	!(SNPInfo$Name \%in\% adjustments$Name) )

#skat
out.skat <- skatMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub)
head(out.skat)

##T1 test
out.t1 <- burdenMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub, mafRange = c(0,0.01))
head(out.t1)

##single snp tests:
out.ss <- singlesnpMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo)
head(out.ss)

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
