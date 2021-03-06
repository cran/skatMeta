\name{skatCohort}
\alias{skatCohort}
\alias{skatFamCohort}
\alias{skatCoxCohort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run SKAT on data from a single cohort.  
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes and organizes the neccesary output to efficiently meta-analyze a linear model SKAT and other tests. Note that the SKAT test is *not* computed by these functions. The output must be passed to one of \code{\link{skatMeta}}, \code{\link{burdenMeta}}, or \code{\link{singlesnpMeta}}.

Unlike the SKAT package which operates on one gene at a time, these functions are intended to operate on many genes, e.g. a whole exome, to facilitate meta analysis of whole genomes or exomes.
}
\usage{
skatCohort(Z, formula, family = gaussian(), SNPInfo=NULL, snpNames = "Name", 
	aggregateBy = "gene", data=parent.frame(), verbose = FALSE)

skatFamCohort(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", 
	data=parent.frame(), fullkins, sparse = TRUE, verbose = FALSE)

skatCoxCohort(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", 
	data=parent.frame(), verbose =FALSE)
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
  \item{snpNames}{
	The field of SNPInfo where the SNP identifiers are found.  Default is 'Name'
	}
  \item{aggregateBy}{
    The field of SNPInfo on which the skat results were aggregated.  Default is 'gene'. For single snps which are intended only for single variant analyses, it is reccomended that they have a unique identifier in this field.  
  }
  \item{data}{
	data frame in which to find variables in the formula
  }	
  \item{verbose}{
  	logical. whether or not to print the progress bar.
  }
  \item{fullkins}{
  	for skatFamCohort: the kinship matrix. See lmekin and the kinship2 package for more details
  }
  \item{sparse}{
  	for skatFamCohort: whether or not to use a sparse Matrix approximation for dense kinship matrices (defaults to TRUE)
  }
}
\details{
	This function computes the neccesary information to meta analyze SKAT analyses: the individual SNP scores, their MAF, and a covariance matrix for each unit of aggregation. Note that the SKAT test is *not* calculated by this function. The output must be passed to one of \code{\link{skatMeta}},\code{\link{burdenMeta}}, or \code{\link{singlesnpMeta}}.
	
A crucial component of SKAT and other region-based tests is a common unit of aggregation accross cohorts. This is given in the SNP information file (argument \code{SNPInfo}), which pairs SNPs to a unit of aggregation (typically a gene). The additional arguments \code{Name} and \code{aggregateBy} specify the columns of the SNP information file which contain these pairings. Note that the column names of the genotype matrix \code{Z} must match the names given in the \code{Name} field.
	
Using \code{skatCohort}, users are strongly recommended to use all SNPs, even if they are monomorphic in your study. This is for two reasons; firstly, monomorphic SNPs provide information about MAF across all studies; without providing the information we are unable to tell if a missing SNP data was monomorphic in a cohort, or simply failed to genotype adequately in that cohort.  Second, even if some SNPs will be filtered out of a particular meta-analysis (e.g., because they are intronic or common) construct ing skatCohort objects describing all SNPs will reduce the workload for subsequent follow-up analyses.


Note: to view results for a single cohort, one can pass a single skatCohort object to a function for meta-analysis.		
}
\value{
	an object of class 'skatCohort'. This is a list, not meant for human consumption, but to be fed to \code{skatMeta()} or another function.
	The names of the list correspond to gene names. Each element in the list contains
	\item{scores}{
		The scores 	(y-yhat)^t g
	}
	\item{cov}{
		The variance of the scores. When no covariates are used, this is the LD matrix.
	}
	\item{n}{
		The number of subjects
	}
	\item{maf}{
		The minor allele frequency
	}
	\item{sey}{
		The residual standard error.
	}


}
\references{
Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics.
		
Chen H, Meigs JB, Dupuis J. Sequence Kernel Association Test for Quantitative Traits in Family Samples.  Genetic Epidemiology. (To appear)

Lin, DY and Zeng, D. On the relative efficiency of using summary statistics versus individual-level data in meta-analysis. Biometrika. 2010.
}
\author{
Arie Voorman, Jennifer Brody
}
\note{
	For \code{skatCoxCohort}, the signed likelihood ratio statistic is used instead of the score, as the score test is anti-conservative for proportional hazards regression. The code for this routine is based on the \code{coxph.fit} function from the \code{survival} package.

Please see the package vignette for more details.
}

\seealso{
\code{\link{skatMeta}}
\code{\link{burdenMeta}}
\code{\link{singlesnpMeta}}
\code{\link{skatOMeta}}
\code{\link{coxph}}
}
\examples{
	
###load example data for two cohorts:
### see ?skatExample	
data(skatExample)

####run on each cohort:
cohort1 <- skatCohort(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
cohort2 <- skatFamCohort(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
	fullkins=kins, data=pheno2)

#### combine results:
##skat
out <- skatMeta(cohort1, cohort2, SNPInfo = SNPInfo)
head(out)

##T1 test
out.t1 <- burdenMeta(cohort1,cohort2, SNPInfo = SNPInfo, mafRange = c(0,0.01))
head(out.t1)

##single snp tests:
out.ss <- singlesnpMeta(cohort1,cohort2, SNPInfo = SNPInfo)
head(out.ss)

\dontrun{
########################
####binary data

cohort1 <- skatCohort(Z=Z1, ybin~1, family=binomial(), SNPInfo = SNPInfo, data =pheno1)
out <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out)


####################
####survival data
cohort1 <- skatCoxCohort(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo = SNPInfo, 
	data =pheno1)
out <- skatMeta(cohort1, SNPInfo = SNPInfo)
head(out)
}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
