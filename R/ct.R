###################################################################################### Read Expression suite file
"readCT"
#' Read files form  Expression suite software output
#' @description Read files form  Expression suite software output
#' @param file Name of the file to read including location
#' @param skip Number of line to skip form reading. Default is 14.
#' @param ... additional parameter to input.
#'
#' @return Data frame. Need to develop in Expression class set.
#'
#' @details The input file can be in txt, csv or xls formate.
#' @author Mohammad Tanvir Ahamed (mashranga@yahoo.com)
#'
#' @examples
#' sfile  <- system.file("extdata", "Ahlam-biopsies_20170920_131310.txt", package = "AsaLab")
#' mdata  <- readCT(file = sfile, skip = 14, header = FALSE)
#'
#' @import tools utils
#'
#'
#' @export
readCT <- function(file, skip = 14,...)
{
  ### Reading file in formate txt, xls or csv
  if(file_ext(file)%in%c("txt","xls","csv"))
  {
    if(file_ext(file)=="txt")
    { targets    <- read.csv(file = file, stringsAsFactors = FALSE, skip=skip, sep = "\t") }
    if(file_ext(file)=="xls")
    { targets    <- gdata::read.xls(xls  = file, sheet=1, skip = skip ) }
    if(file_ext(file)=="csv")
    { targets    <- read.csv(file = file, stringsAsFactors = FALSE, skip = skip ) }
  }
  else
  { stop("Input file should in .txt, .xls or .csv formated", call. = FALSE) }
  targets
}




###################################################################################### Filter expression suite file based on CT value condition
"filterCT"
#' Filter CT value on different condition.
#' @description Filter CT value on different condition. See details.
#' @param file Object output of \link[AsaLab]{readCT}.
#' @param sample Sample column name. Default "Sample.Name"
#' @param target Target column name. Default "Target.Name"
#' @param ct CT column name. Default "Ct"
#' @param ctlimit Limit to keep CT value. Default min= 15 and max = 40.
#' @param undet Value to replace undetermine. Default is NULL. Input value should be numeric.
#' @param ampscore Amplification score to reject the sample. Default 0.6
#' @param phenotype group classification for sample name.
#' @param ctd Absulute difference between ct values to accept. Default is 2. See user documentation for details.
#' @param omit If "TRUE" , sample with TRUE value in omit column will delete. Default is NULL.
#' @param out Name of output file to save the result in CSV formate. Default name "filterCT".
#'
#' @details Details of criteria to delete CT values.
#' For phenotype, The 1st column should be sample name and the column name will be same as the column name of sample in input file.
#' And the second column will be the classification for each sample.The unclassified sample will show NA. If the phonotype
#' dataframe is supplied, in the output dataframe a phenotype colun will show otherwise not.
#'
#' @return This function will return a dataframe of filtered CT values. Also this result will save in CSV formate in current working
#' directory and default name of the file is  filterCT.csv
#'
#' @import psych
#'
#' @examples
#' ###### Example 1
#' # Read expression data
#' sfile  <- "D:\\R Working\\Workflow\\Data\\Asa\\celiac3_20161222_123248_Results_Export.txt"
#' mdata  <- readCT(file = sfile)
#' sname  <- unlist(lapply(strsplit(mdata$Sample.Name," "),"[",1))
#' mdata$Sample.Name <- sname
#' # Read phynotype data
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\Asa\\role-multiple.txt", sep="\t")
#' # Filter CT valued data and add phynotype
#' res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
#'
#'
#' ###### Example 2
#' # Read expression data
#' sfile  <- "D:\\R Working\\Workflow\\Data\\simon_20161102_134801_Results_Export.txt"
#' mdata  <- readCT(file = sfile)
#' # Read phynotype data
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\pheno.txt",sep="\t")
#' # Filter CT valued data and add phynotype
#' res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
#'
#'
#' ###### Example 3
#' # Read expression data
#' sfile  <- "D:\\R Working\\Workflow\\Data\\CT.txt"
#' mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
#' # Read phynotype data
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\phenotypes.txt",sep="\t")
#' # Filter CT valued data and add phynotype
#' res    <- filterCT (file = mdata,sample = "Sample", target = "Detector", phenotype = phen)
#'
#'
#' ####### Example
#' # Read expression data
#' sfile  <- system.file("extdata", "Ahlam-biopsies_20170920_131310.txt", package = "AsaLab")
#' mdata  <- readCT(file = sfile, skip = 14, header = FALSE)
#'
#' # Read phynotype data
#' pfile <- system.file("extdata", "Ahlam-biopsies_20170920_131310_pheno.txt", package = "AsaLab")
#' phen   <- read.csv (file = pfile, sep="\t")
#'
#' # Filter CT valued data and add phynotype
#' res    <- filterCT (file = mdata,sample = "Sample.Name", target = "Target.Name", phenotype = phen)
#'
#' @export
filterCT <- function(file, sample = "Sample.Name",
                     target = "Target.Name", ct = "Ct",
                     ctlimit = c(15,40), ampscore = 0.6,
                     undet = NULL, phenotype = NULL,ctd = 2,
                     omit = NULL, out = "filterCT")
{
  # Delete sample Omit = TRUE
  if(length(omit) > 0 ) file <- file[-which(file$Omit==omit),]
  if(length(omit) ==0 ) file <- file
  # Replace all undetermine CT value with NA
  if(length(undet) > 0)
  {
    if(is.numeric(undet)==FALSE) stop(" The value to replace undetermind should be numeric", call. = TRUE)
    file[,ct][file[,ct]=="Undetermined"]<- undet
  } else
    { file[,ct][file[,ct]=="Undetermined"]<-NA}


  # Replace all CT value with NA outside of limit (ctlimit)
  class(file[,ct])<- "numeric"
  file[,ct][which( file[,ct] < ctlimit[1] | file[,ct] > ctlimit[2])] <- NA
  # Delete sample with amplification score < thrashhold
  #file  <- file[which(file$Amp.Score >= ampscore),]
  # Split data by sample name and target name
  file  <- split(file, file[,sample])
  file  <- lapply(file, function(samplename) { res <- split(samplename, samplename[,target])})
  # Delete the Target name whose all CT == NA
  file1 <- unlist(file, recursive = FALSE)
  file2 <- lapply(file1, function(sam) {res<- all(is.na(sam[,ct]))})
  file  <- file1[which(as.vector(unlist(file2))==FALSE)]
  # Delete sample based on CT value
  myfun_1 <- function(sam)
  {
    sam <- sam[order(sam[,ct]),]
    a   <-as.numeric(sam[,ct])
    # If any vector has missing value
    if ( any(is.na(a))==TRUE )
    {
      # n = 1 non missing
      if (length(a[!is.na(a)]) == 1)
      {
        # Missing value is replace by the single value
        sam[,ct][which(is.na(a))]<- mean(a,na.rm = TRUE)
      }
      # n = 2 non missing
      if (length(a[!is.na(a)]) == 2)
      {
        # Consequtive difference (d) is >  2, delete the sample or if d <= 2, replace missing value by mean
        a1 <- a[!is.na(a)]
        if(abs(as.numeric(dist(a1))) > ctd) sam <-NULL
        else
        {
          sam[,ct][which(is.na(a))]<- mean(a,na.rm = TRUE)
        }
      }
      # n = 3 non missing
      if (length(a[!is.na(a)]) == 3)
      {
        # Difference between all consecutive sample is greater than 2, delete it
        a1 <- a[!is.na(a)]
        d  <- abs(diff(a1))
        if(all(d > ctd))    sam <- NULL
        else
        {
          # d2 > 2d1 or d1 > 2d2
          if(d[2] > 2*d[1])
          {
            sam  <- sam[-which(sam[,ct]==a1[3]),]
            sam[,ct][which(is.na(sam[,ct]))] <- mean(sam[,ct],na.rm = TRUE)
          }
          if(d[1] > 2*d[2])
          {
            sam  <- sam[-which(sam[,ct]==a1[1]),]
            sam[,ct][which(is.na(sam[,ct]))] <- mean(sam[,ct],na.rm = TRUE)
          }
        }
      }
      # n >= 4 non missing
      if (length(a[!is.na(a)]) >= 4)
      {
        # Delete those sample , Ct > median value
        sam <- sam[-which(abs(a-median(a)) > ctd),]
        sam[,ct][which(is.na(sam[,ct]))] <- mean(sam[,ct],na.rm = TRUE)
      }
    }

    # No missing value
    if ( any(is.na(a))==FALSE )
    {
      # n = 2
      if(length(a)==2)
      {
        # Consequtive difference is <= 2
        if(abs(as.numeric(dist(a))) > ctd) sam <-NULL
      }
      # n = 3
      if(length(a) == 3)
      {
        # Difference between all consecutive sample is greater than 2, delete it
        d <- abs(diff(a))
        if(all(d > ctd))    sam <- NULL
        else
        {
          # d2 > 2d1 or d1 > 2d2
          if(d[2] > 2*d[1]) sam <- sam[-3,]
          if(d[1] > 2*d[2]) sam <- sam[-1,]
        }
      }
      # n = 4
      if(length(a) >= 4)
      {
        # Take madian of 4 and delete which is greater than median
        sam <- sam[which(abs(a-median(a)) < ctd),]
      }
    }
    sam
  }

  file <- lapply(file, myfun_1)

  # Get mean on ct value on every splied section by sample and target
  mufun_2 <- function(i)
  {
    res <- data.frame(t(c(unique(i[,sample]), unique(i[,target]), round(mean(as.numeric(i[,ct])),3 ))))
    if (length(res) >= 2 ) names(res) <- c(sample, target , "Ct.Mean" )
    else res <- NULL
    res
  }
  file <- lapply(file,mufun_2)
  file <- do.call(rbind, file)

  ###### Reshape data
  file <- reshape(file, v.names = "Ct.Mean", idvar = sample,  timevar = target, direction = "wide")
  nam   <- c(sample, unlist(lapply(strsplit(names(file)[-1],"Ct.Mean."),"[[",2)))
  names(file) <- nam
  samname <- file[,sample]
  file[,sample] <- NULL
  file <- apply(file,2, function(i) as.numeric(as.character(i)) )
  rownames(file) <- samname
  file1<- file

  ###### Add phenotype data [Only if phenotype data supplied]

  if (length(phenotype) >  0 )
  {
    if(ncol(phenotype) > 2 ) stop(" Phoneotype file should contain 2 column. Check the file.", call. = FALSE)
    file1 <- data.frame(file1)
    file1 <- merge(phenotype,file1, by.x= colnames(phenotype)[1], by.y = "row.names")
    file1
  }

  if (length(phenotype) ==  0 )
  {
    sample <- rownames(file1)
    file1  <- cbind(sample,file1)
    file1  <- data.frame(file1)
    rownames(file1) <- NULL
  }
  write.csv(x = file1, file = paste0(out,".csv"), row.names = FALSE)

  file1
}



###################################################################################### Analyse CT values
"analyzeCT"
#' Analyze CT value on different condition.
#' @description Analyze CT values on defferent condition
#'
#' @param file data in Matrix. Output of \link[AsaLab]{filterCT} function.
#' @param skip Number of column to skip form 1st column in input file. Default is 2.
#' @param hkgene House keeping gene.
#' @param control Name of phynotype that will assign as control. Default is "CONTROL".
#' @param del Percentage of missing value for which any gene or sample will be excuded form analysis. Default is 0.7 (70 percentage).
#' @param missCT If missing value will replace with 40 or otheres. Default is NULL.
#' @param gylim Define Upper and lower limit for the boxplot. Default is c(-30,30).
#'
#' @return This function will return a list of deleted genes, deleted sample, delta CT values based on supplied housekeeping genes
#' and delta CT values based average of supplied housekeeping genes.
#' CT values will be saved in a csv file after deleteing the gene and sample with preferred percentage of missing values. Also save results
#' of the boxplot of delta CT values and p-values for t-test based on both individual and average of
#' housekeeping genes. Boxplot and p values will save in pdf formate and delta CT values based on both individual and average
#' housekeeping genes will save in CSV formate in current working directory.
#' @import impute
#' @import psych
#' @import lattice
#' @import graphics
#' @import grDevices
#' @import stats
#'
#' @examples
#' ####### Example 1
#' ## Read expression data #Sample : 384
#' sfile  <- "D:\\R Working\\Workflow\\Data\\Asa\\celiac3_20161222_123248_Results_Export.txt"
#' mdata  <- readCT(file = sfile)
#' sname  <- unlist(lapply(strsplit(mdata$Sample.Name," "),"[",1))
#' mdata$Sample.Name <- sname
#' ## Read phynotype data #Sample : 446
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\Asa\\role-multiple.txt", sep="\t")
#' ## Filter CT valued data and add phynotype #Sample : 330
#' res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
#' ## Analyze CT values
#' hkg    <- c("YWHAZ_047","GUSB_627","HPRT1_909")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="CONTROL",del = 0.7, missCT = NULL)
#'
#'
#'
#' ###### Example 2
#' ## Read expression data #Sample : 80
#' sfile  <- "D:\\R Working\\Workflow\\Data\\simon_20161102_134801_Results_Export.txt"
#' mdata  <- readCT(file = sfile)
#' ## Read phynotype data #Sample : 80
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\pheno.txt",sep="\t")
#' ## Filter CT valued data and add phynotype #Sample : 80
#' res    <- filterCT (file = mdata,omit = "True",phenotype = phen)
#' ## Analyze CT values
#' hkg    <- c("YWHAZ","GUSB","HSPB1")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="control",del = 0.7, missCT = 40)
#'
#'
#'
#' ###### Example 3
#' ## Read expression data #Sample : 55
#' sfile  <- "D:\\R Working\\Workflow\\Data\\CT.txt"
#' mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
#' ## Read phynotype data  #Sample : 114
#' phen   <- read.csv (file="D:\\R Working\\Workflow\\Data\\phenotypes.txt",sep="\t")
#' ## Filter CT valued data and add phynotype #Sample : 54
#' res    <- filterCT (file = mdata,sample = "Sample", target = "Detector", phenotype = phen)
#' hkg    <- c("ZMYM2","GUSB_627","HFE2")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="0",del = 0.7, missCT = 40)
#'
#'
#' ###### Example 4
#' ## Read expression data #Sample : 135
#' sfile  <- "Z:\\TANVIR\\George\\Nasal Polys-George_20160808_150714_Results_Export (2).txt"
#' mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
#' ## Read phynotype data #Sample : 72
#' phen   <- read.csv (file="Z:\\TANVIR\\George\\phenotypes.txt" ,sep="\t")
#' ## Filter CT valued data and add phynotype #Sample : 71
#' res    <- filterCT (file = mdata,sample = "Sample.Name", target = "Target.Name", phenotype = phen)
#' hkg    <- c("HPRT1","YWHAZ")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkg,control="0",del = 0.1, missCT = 40)
#'
#'
#' ###### Example 5
#' ## Read expression data #Sample : 135
#' sfile  <- "Z:\\TANVIR\\George\\control sample\\CD for analysis.txt"
#' mdata  <- readCT(file = sfile, skip = 0, header = FALSE)
#' ## Read phynotype data  #Sample : 72
#' phen   <- read.csv (file="Z:\\TANVIR\\George\\control sample\\phenotypes.txt" ,sep="\t")
#' ## Filter CT valued data and add phynotype #Sample : 71
#' res    <- filterCT(file=mdata,sample="Sample.Name",
#'                    target="Target.Name",undet=40,ctd=1,phenotype=phen)
#' hkge    <- c("HPRT1.909","YWHAZ.047")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkge,control="0",del = 0.8, missCT = NULL)
#'
#'
#' ##### Example
#' #' # Read expression data
#' sfile  <- system.file("extdata", "Ahlam-biopsies_20170920_131310.txt", package = "AsaLab")
#' mdata  <- readCT(file = sfile, skip = 14, header = FALSE)
#'
#' # Read phynotype data
#' pfile <- system.file("extdata", "Ahlam-biopsies_20170920_131310_pheno.txt", package = "AsaLab")
#' phen   <- read.csv (file = pfile, sep="\t")
#'
#' # Filter CT valued data and add phynotype
#' res    <- filterCT (file = mdata,sample = "Sample.Name", target = "Target.Name", phenotype = phen)
#'
#' hkge    <- c("DUSP1","HPRT1","YWHAZ")
#' CTres  <- analyzeCT (file = res, skip=2, hkgene = hkge,control="CONTROL",del = 0.9, missCT = NULL)
#'
#'
#'
#' @export
analyzeCT <- function(file, skip = 2, hkgene ,control="CONTROL", del = 0.7, missCT = NULL,gylim = c(-30,30))
{

  if (class(file[,1])=="integer"){ class(file[,1])<-"character" }
  if (class(file[,2])=="integer"){ class(file[,2])<-"character" }

  ###### Check house keepimng gene

  if(any(hkgene%in%colnames(file)==FALSE)) {stop("\nError : House keeping gene not in data. Check house keeping gene name.\n","- Availabel gene name:\n",paste0(colnames(file)[-c(1:2)],sep = "; "),"\n- Supplied Housekeeping gene:\n",paste0(hkgene,sep = "; "),call. = TRUE)}


  ###### Separate CT values and sample information
  sample <- file [,c(1:skip)]
  ct     <- data.frame(as.matrix(file [,-c(1:skip)]))

  ###### Delete gene with 70% missing value and Missing value imputation for house keeping gene
  # Delete Gene for missing value
  if (length(which(colMeans(is.na(file)) > del)) > 0 )
  {
    file    <- ct[, -which(colMeans(is.na(ct)) > del) ]
    #sample <- sample[, -which(colMeans(is.na(ct)) > del) ]
    delGene <- colnames(ct)[which(colMeans(is.na(ct)) > del)]
    message("Deleted Genes names for ",del*100,"% of missing values : \n", paste0(delGene, sep = "; "))
  } else
  {
    file   <- ct
    sample <- sample
    delGene <- NA
    message("Deleted Genes names for ",del*100,"% of missing values : \n", delGene)
  }
  # Delete sample for missing value
  if (length(which(rowMeans(is.na(file)) > del)) > 0 )
  {
    file1      <- file[-which(rowMeans(is.na(file)) > del), ]
    delSample  <- unique(as.vector(sample[,1][which(rowMeans(is.na(file)) > del)])) # Deleted sample
    sample1    <- sample[-which(rowMeans(is.na(file)) > del),]

    file       <- file1
    sample     <- sample1

    message("\nDeleted sample names for ",del*100,"% of missing values : \n", paste0(delSample, sep = "; "))
  } else
  {
    file   <- file
    sample <- sample
    delSample <- NA
    message("\nDeleted sample names for ",del*100,"% of missing values : \n", NA)
  }

  # Impute missing value
  #require(impute)
  #file   <- suppressWarnings(t(impute.knn(t(file))$data))

  # Replace missing value with 40 or not
  if(length(missCT) == 0)
  {
    file <- file
  } else
  {
    file[is.na(file)]<- missCT
    file <- file
  }

  ###### Save file
  write.csv(x = cbind(sample,file), file = "filterCT_1.csv", row.names = FALSE)

  ###### House keeping gene : Delta ct valuses
  # Availabel house keeping gene
  hkg      <- as.list(hkgene[hkgene%in%colnames(file)])
  hkg_ge   <- file[,match(hkg, colnames(file))]

  hkg_file <- data.frame(file[which(rowMeans(is.na(hkg_ge))   <= (length(hkg)-1)/length(hkg)), ])   ## Filter missing values
  hkg_sa   <- data.frame(sample[which(rowMeans(is.na(hkg_ge)) <= (length(hkg)-1)/length(hkg)), ])   ## Filter missing values
  hkg_ge   <- data.frame(hkg_ge[which(rowMeans(is.na(hkg_ge)) <= (length(hkg)-1)/length(hkg)), ])   ## Filter missing values

  #write.csv(x = cbind(sample,file[,match(hkg, colnames(file))]), file = "hkg.csv", row.names = FALSE)
  write.csv(x = cbind(hkg_sa,hkg_ge), file = "hkg.csv", row.names = FALSE)

  # Average CT values form house keeping gene
  #avgHKG <- rowMeans(file[,match(hkg, colnames(file))], na.rm = TRUE)
  #avgCT  <- round(file-avgHKG,3)
  #avgCT  <- data.frame(cbind(sample,avgCT))
  #write.csv(x = avgCT, file = "avg_deltaCT.csv", row.names = FALSE)

  avgHKG <- rowMeans(t(impute.knn(t(hkg_ge))$data), na.rm= TRUE)
  #avgCT  <- round(hkg_file-avgHKG,3)
  avgCT  <- round(avgHKG-hkg_file,3)
  avgCT  <- data.frame(cbind(hkg_sa,avgCT))
  write.csv(x = avgCT, file = "avg_deltaCT.csv", row.names = FALSE)

  # Form individual house keeping gene
  file   <- lapply(hkg,function(i)
  {
    res0 <- file[,!(colnames(file) == i)]
    res1 <- file[,i]
    #res1 <- file[,i,drop=FALSE]
    res2 <- round(res1-res0,3)
    res3 <- data.frame(cbind(sample,res2))
    write.csv(x = res3, file = paste0(i,"_deltaCT",".csv"), row.names = FALSE)
    res3
  }
  )
  file$avCT <- avgCT

  names(file) <- c(hkg, "avg_deltaCT")
  delCT       <- file

  ###### Split data by phinotype
  file <- lapply(file, function(i) split(i, i[,2]))

  ###### Reshape data
  #library(reshape)
  file  <- suppressMessages (lapply(file, reshape::melt))
    file <- lapply(file, function(i)
  {
    res <- split(i, i[,3])
    res
  })
  #print(file[[1]][[1]])


  ###### Ploting Gene
  res <- suppressWarnings( lapply(seq_along(file), function(i)
    {
    ## Name of refgene
    nam <- names(file[i])
    pdf(paste0(nam,".pdf"))

    n <- length(file[[i]])
    for(j in 1:n )
    {
      v <- file[[i]][[j]]
      # Change levels to put CONTROL group at first
      #print(paste0(i, " and ", j))
      v[,2] <- factor(v[,2], levels= c(control,setdiff(unique(v[,2]),control)))
      v     <- v[order(v[,2]), ]

      # Scalling data based on control group
      v1 <- split(v,v[,2])
      cn <- which(names(v1)==control)
      mn <- mean(as.vector(v1[[cn]][,4]),na.rm = TRUE)
      v2 <- v
      v2[,4] <- v2[,4]-mn
      v<-v2

      # Group mean
      po  <- split(v,v[,2])
      gmn <- lapply(po, function(i)
      {
        k <- mean(as.vector(i[,4]),na.rm = TRUE)
      })
      gmn <- unlist(gmn)

      ####### Plot bar plot and p-value
      op <- par()
      par(mfrow = c(1, 2))

      # P-value
      g<- v
      colnames(g)[2]<- "CD"
      kp<- combn(as.character(unique(g$CD)), 2,
                 FUN = function(y)
                 {
                   #print(y)
                   gk <- droplevels(g[g$CD %in% y ,])
                   if(any(rowSums(table(gk$CD,gk$value))<=1))
                    {
                     res <- list(0)
                     res$p.value <- NA
                     res
                    } else{
                   res<- t.test(value ~ CD, g[g$CD %in% y ,])
                   res
                    }
                 }
                 , simplify = FALSE)

      pp<- unlist(lapply(kp, function(i) i$p.value))

      combs  <- combn(as.character(unique(g$CD)), 2)
      N      <- pp
      inputs <- as.character(unique(g$CD))
      out2 <- N
      class(out2) <- "dist"
      attr(out2, "Labels") <- as.character(inputs)
      attr(out2, "Size") <- length(inputs)
      attr(out2, "Diag") <- attr(out2, "Upper") <- FALSE
      out2 <- round(as.matrix(out2),3)

      out2[upper.tri(out2,diag = TRUE)] <- "-"

      if (ncol(out2) > 2)
      {
      out2 <- out2[,colSums(is.na(out2)) < nrow(out2)]
      out2 <- out2[rowSums(is.na(out2)) < ncol(out2),]
      }

      #require(lattice)
      myPanel <- function(x, y, z, ...) {
        panel.levelplot(x,y,z,...)
        panel.text(x, y,  out2[cbind(x,y)],cex = 0.3) ## use handy matrix indexing
      }

      p1<- suppressWarnings( levelplot(out2,panel =myPanel, scales=list(x=list(rot=90)),
                     xlab ="", ylab ="", main = "p-value",sub = "NA = Either 1 or 0 observation\n[No t-test has performed] ",
                     at = seq(0,0.05,by =  0.001),
                     col.regions = c(heat.colors(51))))


      plot(1,axes=FALSE, col = "white", xlab = "", ylab="")
      print(p1, position = c(0, 0, 0.5, 1), more = TRUE)


      # Bar plot
      if(length(names(table(v[,2]))) ==2 )
      {
        #boxplot(v[which(v[,2]==names(table(v[,2]))[1]),][,4], xlim = c(0.5, 2+0.5),ylim = range (v[,4],na.rm = TRUE),col="white",
        #        boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), ylab = "Fold",
        #        main = paste0("Ref gene: ", nam," \nGene: ", names(file[[i]][j] ))) #invisible boxes

        boxplot(v[which(v[,2]==names(table(v[,2]))[1]),][,4], xlim = c(0.5, 2+0.5),ylim = c(gylim[1],gylim[2]),col="white",
                boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1), ylab = "Fold",
                main = paste0("Ref gene: ", nam," \nGene: ", names(file[[i]][j] ))) #invisible boxes

        boxplot(v[which(v[,2]==names(table(v[,2]))[1]),][,4], xaxt = "n", add = TRUE,at = 1) #shift these left by -0.15
        points(rep(1, length(v[which(v[,2]==names(table(v[,2]))[1]),][,4] )),v[which(v[,2]==names(table(v[,2]))[1]),][,4],type = "p",cex=0.5, col = "orange")
        points(1,gmn[1],col="blue",lwd = 2)

        boxplot(v[which(v[,2]==names(table(v[,2]))[2]),][,4], xaxt = "n", add = TRUE,at = 2) #shift these left by -0.15
        points(rep(2,length(v[which(v[,2]==names(table(v[,2]))[2]),][,4]) ),v[which(v[,2]==names(table(v[,2]))[2]),][,4],type = "p",cex=0.5, col = "orange")
        points(2,gmn[2],col="blue",lwd = 2)

        axis(1, at=c(1,2), labels = names(table(v[,2])) )


      }  else
        {
          #suppressWarnings(plot(v[,2],v[,4], main = paste0("Ref gene: ", nam," \nGene: ", names(file[[i]][j] )), ylab="Fold", xlab=NA,las = 2, cex.axis = 0.7,
          # ylim = c( ifelse(is.infinite(min(v[,4],na.rm= TRUE))==TRUE,-10,min(v[,4],na.rm= TRUE))  ,
          #           ifelse(is.infinite(max(v[,4],na.rm= TRUE))==TRUE,10,max(v[,4],na.rm= TRUE))) ))

          suppressWarnings(plot(v[,2],v[,4], main = paste0("Ref gene: ", nam," \nGene: ", names(file[[i]][j] )), ylab="Fold", xlab=NA,las = 2, cex.axis = 0.7,
                                ylim = c(gylim[1],gylim[2]) ))



      points(v[,2],v[,4],type = "p",cex=0.5, col = "orange")
      points(gmn, col = "blue", lwd = 2)
      }

      suppressWarnings (par(op))
      }

    dev.off()

    }))

  res <- list (del_gene = delGene,del_sample = delSample, delta_CT = delCT)
}


