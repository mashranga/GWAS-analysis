## ---- echo=FALSE---------------------------------------------------------
library(AsaLab)
data(ex_data_1)
head(ex_data_1)

## ---- echo=FALSE---------------------------------------------------------
library(AsaLab)
data(ex_phenotype_1)
head(ex_phenotype_1)

## ------------------------------------------------------------------------
## CT values form ExpressionSuite output
data(ex_data_1)
mdata <- ex_data_1
## Phenotype file
data(ex_phenotype_1)
phen <- ex_phenotype_1

## Apply CT value filtering criteria
# Without phenotype file input
res  <- filterCT(file = mdata, phenotype = NULL, out = "filterCT")
head(res)
# With phenotype file input
res  <- filterCT(file = mdata, phenotype = phen, out = "filterCT")
head(res)

## ---- message=FALSE------------------------------------------------------
## CT values form ExpressionSuite output
data(ex_data_1)
mdata <- ex_data_1
## Phenotype file
data(ex_phenotype_1)
phen <- ex_phenotype_1
## Filter CT values 
res  <- filterCT(file = mdata, phenotype = phen, out = "filterCT")
## House keeping gene
hkg    <- c("YWHAZ_047","GUSB_627","HPRT1_909")
CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="control",del = 0.7, missCT = NULL)

## ---- echo=FALSE---------------------------------------------------------
p1 <- CTres$del_gene
p1

## ---- echo=FALSE---------------------------------------------------------
p2 <- CTres$del_sample
p2

## ---- echo=FALSE---------------------------------------------------------
p3<- lapply(CTres$delta_CT,head)
p3

