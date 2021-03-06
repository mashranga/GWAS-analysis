% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ct.R
\name{analyzeCT}
\alias{analyzeCT}
\title{Analyze CT value on different condition.}
\usage{
analyzeCT(file, skip = 2, hkgene, control = "CONTROL", del = 0.7,
  missCT = NULL, gylim = c(-30, 30))
}
\arguments{
\item{file}{data in Matrix. Output of \link[AsaLab]{filterCT} function.}

\item{skip}{Number of column to skip form 1st column in input file. Default is 2.}

\item{hkgene}{House keeping gene.}

\item{control}{Name of phynotype that will assign as control. Default is "CONTROL".}

\item{del}{Percentage of missing value for which any gene or sample will be excuded form analysis. Default is 0.7 (70 percentage).}

\item{missCT}{If missing value will replace with 40 or otheres. Default is NULL.}

\item{gylim}{Define Upper and lower limit for the boxplot. Default is c(-30,30).}
}
\value{
This function will return a list of deleted genes, deleted sample, delta CT values based on supplied housekeeping genes
and delta CT values based average of supplied housekeeping genes.
CT values will be saved in a csv file after deleteing the gene and sample with preferred percentage of missing values. Also save results
of the boxplot of delta CT values and p-values for t-test based on both individual and average of
housekeeping genes. Boxplot and p values will save in pdf formate and delta CT values based on both individual and average
housekeeping genes will save in CSV formate in current working directory.
}
\description{
Analyze CT values on defferent condition
}
\examples{
####### Example 1
## Read expression data #Sample : 384
sfile  <- "D:\\\\R Working\\\\Workflow\\\\Data\\\\Asa\\\\celiac3_20161222_123248_Results_Export.txt"
mdata  <- readCT(file = sfile)
sname  <- unlist(lapply(strsplit(mdata$Sample.Name," "),"[",1))
mdata$Sample.Name <- sname
## Read phynotype data #Sample : 446
phen   <- read.csv (file="D:\\\\R Working\\\\Workflow\\\\Data\\\\Asa\\\\role-multiple.txt", sep="\\t")
## Filter CT valued data and add phynotype #Sample : 330
res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
## Analyze CT values
hkg    <- c("YWHAZ_047","GUSB_627","HPRT1_909")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="CONTROL",del = 0.7, missCT = NULL)



###### Example 2
## Read expression data #Sample : 80
sfile  <- "D:\\\\R Working\\\\Workflow\\\\Data\\\\simon_20161102_134801_Results_Export.txt"
mdata  <- readCT(file = sfile)
## Read phynotype data #Sample : 80
phen   <- read.csv (file="D:\\\\R Working\\\\Workflow\\\\Data\\\\pheno.txt",sep="\\t")
## Filter CT valued data and add phynotype #Sample : 80
res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
## Analyze CT values
hkg    <- c("YWHAZ","GUSB","HSPB1")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="control",del = 0.7, missCT = 40)



###### Example 3
## Read expression data #Sample : 55
sfile  <- "D:\\\\R Working\\\\Workflow\\\\Data\\\\CT.txt"
mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
## Read phynotype data  #Sample : 114
phen   <- read.csv (file="D:\\\\R Working\\\\Workflow\\\\Data\\\\phenotypes.txt",sep="\\t")
## Filter CT valued data and add phynotype #Sample : 54
res    <- filterCT (file = mdata,sample = "Sample", target = "Detector", phenotype = phen)
hkg    <- c("ZMYM2","GUSB_627","HFE2")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="0",del = 0.7, missCT = 40)


###### Example 4
## Read expression data #Sample : 135
sfile  <- "Z:\\\\TANVIR\\\\George\\\\Nasal Polys-George_20160808_150714_Results_Export (2).txt"
mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
## Read phynotype data #Sample : 72
phen   <- read.csv (file="Z:\\\\TANVIR\\\\George\\\\phenotypes.txt" ,sep="\\t")
## Filter CT valued data and add phynotype #Sample : 71
res    <- filterCT (file = mdata,sample = "Sample.Name", target = "Target.Name", phenotype = phen)
hkg    <- c("HPRT1","YWHAZ")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="0",del = 0.1, missCT = 40)


###### Example 5
## Read expression data #Sample : 135
sfile  <- "Z:\\\\TANVIR\\\\George\\\\control sample\\\\CD for analysis.txt"
mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
## Read phynotype data  #Sample : 72
phen   <- read.csv (file="Z:\\\\TANVIR\\\\George\\\\control sample\\\\phenotypes.txt" ,sep="\\t")
## Filter CT valued data and add phynotype #Sample : 71
res    <- filterCT(file=mdata,sample="Sample.Name",
                   target="Target.Name",undet=40,ctd=1,phenotype=phen)
hkge    <- c("HPRT1.909","YWHAZ.047")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkge,control="0",del = 0.8, missCT = NULL)


##### Example
#' # Read expression data
sfile  <- system.file("extdata", "Ahlam-biopsies_20170920_131310.txt", package = "AsaLab")
mdata  <- readCT(file = sfile, skip = 14, header = FALSE)

# Read phynotype data
pfile <- system.file("extdata", "Ahlam-biopsies_20170920_131310_pheno.txt", package = "AsaLab")
phen   <- read.csv (file = pfile, sep="\\t")

# Filter CT valued data and add phynotype
res    <- filterCT (file = mdata,sample = "Sample.Name", target = "Target.Name", phenotype = phen)

hkge    <- c("DUSP1","HPRT1","YWHAZ")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkge,control="CONTROL",del = 0.9, missCT = NULL)



}
