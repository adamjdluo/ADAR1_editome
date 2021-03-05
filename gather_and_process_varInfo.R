## The script is designed for further variant processing, integrating information and generate report 
##
## Environment
dp_cut <- 5 # depth >=5
ac_cut <- 2 # alternative allele count >=2
num_cut <- 3 # in at least three sample
load("hg19_repeatMask.RData")
load("path to variants RData")
#
base_dir <- "filePath to the project directory"
res_dir <- file.path(base_dir,"result directory")
dir.create(res_dir)
##
## Process variant annotation
annoTXT <- file.path(base_dir,"annoTab.txt")
annoTab <- read.table(annoTXT,header = FALSE,sep="\t",stringsAsFactors = FALSE)
names(annoTab) <- c("variant","geneSymbol","SNP138","geneFunc","exonFunc","AAChange","Description")
annoTab <- dplyr::select(annoTab,variant,geneSymbol,SNP138,geneFunc,exonFunc,AAChange,Description)
annoTab$geneFunc[grepl(annoTab$geneFunc,pattern = "ncRNA")] <- "ncRNA"
annoTab$geneFunc[grepl(annoTab$geneFunc,pattern = "downstream")|
                   grepl(annoTab$geneFunc,pattern = "upstream")] <- "intergenic"
annoTab$geneFunc[grepl(annoTab$geneFunc,pattern = "UTR5;UTR3")] <- "UTR3"
annoTab$geneFunc[grepl(annoTab$geneFunc,pattern = "exonic;splicing")] <- "exonic"
annoTab$geneFunc[grepl(annoTab$geneFunc,pattern = "exonic")] <- "CDS"
##
## Integrate information into single table
# depth of each sample
varDP <- data.frame(variant=rownames(matDP))
varDP <- cbind(varDP,matDP)
names(varDP) <- c("variant",paste0("DP_",names(varDP)[-1]))
varDP$variant <- as.vector(varDP$variant)
#
# alleleic frequency (AF) for each sample
varAF <- data.frame(variant=rownames(matAF))
varAF <- cbind(varAF,matAF)
names(varAF) <- c("variant",paste0("AF_",names(varAF)[-1]))
varAF$variant <- as.vector(varAF$variant)
#
# alleleic counts (AC) for each samplle
varAC <- data.frame(variant=rownames(matAC))
varAC <- cbind(varAC,matAC)
names(varAC) <- c("variant",paste0("AC_",names(varAC)[-1]))
varAC$variant <- as.vector(varAC$variant)
#
# Merge DP, AF, and AC into single table
bigTab <- merge(varDP,varAF,id='variant')
bigTab <- merge(bigTab,varAC,id='variant')
#
# other information for variants
rd_vcf_Alu <- rd[unique(queryHits(GenomicRanges::findOverlaps(rd,rd_Alu)))] # list of Alu elements
var_inAlu <- names(rd_vcf_Alu)
bigTab$inAlu <- as.integer(bigTab$variant %in% var_inAlu) # decide the given variants in Alu elements or not
bigTab$chr <- gsub("(.*):(.*)_(.*)","\\1",bigTab$variant) # extract chromosome
bigTab$position <- gsub("(.*):(.*)_(.*)","\\2",bigTab$variant) # extract position
bigTab$ntCh <- gsub("(.*):(.*)_(.*)","\\3",bigTab$variant) # nucleotide changes
#
# separate samples with/without IFN treatment
bigTab_ori <- bigTab[,c(names(bigTab)[!grepl(names(bigTab),pattern = "-2")])] # remove IFN samples
bigTab_IFN <- bigTab[,c('variant',names(bigTab)[grepl(names(bigTab),pattern = "-2")],'inAlu','ntCh')] # retend IFN sample only
#
bigTab_ori <- merge(bigTab_ori,annoTab,id='variant',all.x=TRUE)
bigTab_IFN <- merge(bigTab_IFN,annoTab,id='variant',all.x=TRUE)
##
## Select candidates in non-IFN
tab_pf <- bigTab_ori
tab_pf$dpCut_C <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_C")]]>=dp_cut)
tab_pf$dpCut_F <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_F")]]>=dp_cut)
tab_pf$dpCut_X <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_X")]]>=dp_cut)
tab_pf$dpCut_Z <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_Z")]]>=dp_cut)
#
tab_pf$acCut_C <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_C")]]>=ac_cut)
tab_pf$acCut_F <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_F")]]>=ac_cut)
tab_pf$acCut_X <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_X")]]>=ac_cut)
tab_pf$acCut_Z <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_Z")]]>=ac_cut)
#
tab_pf$psCut_C <- as.integer(tab_pf$dpCut_C==3&tab_pf$acCut_C==3)
tab_pf$psCut_F <- as.integer(tab_pf$dpCut_F==3&tab_pf$acCut_F==3)
tab_pf$psCut_X <- as.integer(tab_pf$dpCut_X==3&tab_pf$acCut_X>=1) # either one pass the
tab_pf$psCut_Z <- as.integer(tab_pf$dpCut_Z==3&tab_pf$acCut_Z==3)
#
tab_pf$psCut_CX <- as.integer(tab_pf$psCut_C==1&tab_pf$psCut_X==0)
tab_pf$psCut_FX <- as.integer(tab_pf$psCut_F==1&tab_pf$psCut_X==0)
tab_pf$psCut_ZX <- as.integer(tab_pf$psCut_Z==1&tab_pf$psCut_X==0)
#
tab_filt_ct <- tab_pf[tab_pf$psCut_CX==1|tab_pf$psCut_FX==1|tab_pf$psCut_ZX==1,]
#
# filter by AF changes
cL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_C")]
fL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_F")]
xL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_X")]
zL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_Z")]
#
for(k in 1:length(tab_filt_ct$variant)){
  cnum <- unlist(tab_filt_ct[k,cL])
  fnum <- unlist(tab_filt_ct[k,fL])
  xnum <- unlist(tab_filt_ct[k,xL])
  znum <- unlist(tab_filt_ct[k,zL])
  #
  if(sum(xnum)==0){cutX <- 0}else{cutX<-mean(xnum)+2*sd(xnum)}
  #
  tab_filt_ct$sigCX[k] <- as.integer(sum(as.integer(cnum>cutX))==3)
  tab_filt_ct$sigFX[k] <- as.integer(sum(as.integer(fnum>cutX))==3)
  tab_filt_ct$sigZX[k] <- as.integer(sum(as.integer(znum>cutX))==3)}
#
# grouping variants
tab_filt_ct_AF <- tab_filt_ct[tab_filt_ct$sigCX==1|tab_filt_ct$sigFX==1|tab_filt_ct$sigZX==1,]
#
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 4
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 5
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 6
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 7
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 8
#
tab_filt_ct_AF_inAlu <- rbind(tab_filt_ct_AF[tab_filt_ct_AF$inAlu==1&tab_filt_ct_AF$ntCh=="A/G",], # select A/G mutations and the position in Alu element
                              tab_filt_ct_AF[tab_filt_ct_AF$inAlu==1&tab_filt_ct_AF$ntCh=="T/C",]) # select T/C mutations and the position in Alu element
tab_filt_ct_AF_inAlu <- tab_filt_ct_AF_inAlu[!grepl(tab_filt_ct_AF_inAlu$SNP138,pattern = 'rs'),] #select variants not in SNP138 database
#
rio::export(list("filtByCount"=tab_filt_ct,"filtByAF"=tab_filt_ct_AF,"inAlu"=tab_filt_ct_AF_inAlu),
            file = file.path(res_dir,"var_cand_nIFN.xlsx"))
##
## Select candidates in IFN
tab_pf <- bigTab_IFN
tab_pf$dpCut_C <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_C")]]>=dp_cut)
tab_pf$dpCut_F <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_F")]]>=dp_cut)
tab_pf$dpCut_X <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_X")]]>=dp_cut)
tab_pf$dpCut_Z <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "DP_Z")]]>=dp_cut)
#
tab_pf$acCut_C <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_C")]]>=ac_cut)
tab_pf$acCut_F <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_F")]]>=ac_cut)
tab_pf$acCut_X <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_X")]]>=ac_cut)
tab_pf$acCut_Z <- rowSums(tab_pf[,names(tab_pf)[grepl(names(tab_pf),pattern = "AC_Z")]]>=ac_cut)
#
tab_pf$psCut_C <- as.integer(tab_pf$dpCut_C==3&tab_pf$acCut_C==3)
tab_pf$psCut_F <- as.integer(tab_pf$dpCut_F==3&tab_pf$acCut_F==3)
tab_pf$psCut_X <- as.integer(tab_pf$dpCut_X==3&tab_pf$acCut_X>=1) # either one pass the
tab_pf$psCut_Z <- as.integer(tab_pf$dpCut_Z==3&tab_pf$acCut_Z==3)
#
tab_pf$psCut_CX <- as.integer(tab_pf$psCut_C==1&tab_pf$psCut_X==0)
tab_pf$psCut_FX <- as.integer(tab_pf$psCut_F==1&tab_pf$psCut_X==0)
tab_pf$psCut_ZX <- as.integer(tab_pf$psCut_Z==1&tab_pf$psCut_X==0)
#
tab_filt_ct <- tab_pf[tab_pf$psCut_CX==1|tab_pf$psCut_FX==1|tab_pf$psCut_ZX==1,]
#
# filter by AF changes
cL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_C")]
fL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_F")]
xL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_X")]
zL <- names(tab_filt_ct)[grepl(names(tab_filt_ct),pattern = "AF_Z")]
#
for(k in 1:length(tab_filt_ct$variant)){
  cnum <- unlist(tab_filt_ct[k,cL])
  fnum <- unlist(tab_filt_ct[k,fL])
  xnum <- unlist(tab_filt_ct[k,xL])
  znum <- unlist(tab_filt_ct[k,zL])
  #
  if(sum(xnum)==0){cutX <- 0}else{cutX<-mean(xnum)+2*sd(xnum)}
  #
  tab_filt_ct$sigCX[k] <- as.integer(sum(as.integer(cnum>cutX))==3)
  tab_filt_ct$sigFX[k] <- as.integer(sum(as.integer(fnum>cutX))==3)
  tab_filt_ct$sigZX[k] <- as.integer(sum(as.integer(znum>cutX))==3)}
#
# grouping variants
tab_filt_ct_AF <- tab_filt_ct[tab_filt_ct$sigCX==1|tab_filt_ct$sigFX==1|tab_filt_ct$sigZX==1,]
#
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 1
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 2
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 3
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 4
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 5
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&tab_filt_ct_AF$sigZX==1] <- 6
tab_filt_ct_AF$Grp[tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 7
tab_filt_ct_AF$Grp[!tab_filt_ct_AF$sigCX==1&!tab_filt_ct_AF$sigFX==1&!tab_filt_ct_AF$sigZX==1] <- 8
#
tab_filt_ct_AF_inAlu <- rbind(tab_filt_ct_AF[tab_filt_ct_AF$inAlu==1&tab_filt_ct_AF$ntCh=="A/G",], # select A/G mutations and the position in Alu element
                              tab_filt_ct_AF[tab_filt_ct_AF$inAlu==1&tab_filt_ct_AF$ntCh=="T/C",]) # select T/C mutations and the position in Alu element
tab_filt_ct_AF_inAlu <- tab_filt_ct_AF_inAlu[!grepl(tab_filt_ct_AF_inAlu$SNP138,pattern = 'rs'),] #select variants not in SNP138 database
#
rio::export(list("filtByCount"=tab_filt_ct,"filtByAF"=tab_filt_ct_AF,"inAlu"=tab_filt_ct_AF_inAlu),
            file = file.path(res_dir,"var_cand_IFN.xlsx"))
##
## Write report
file.copy(file.path(base_dir,"report_rmd_v2.Rmd"),res_dir,recursive = TRUE)
file_ori <- file.path(res_dir,"var_cand_nIFN.xlsx")
rmarkdown::render(file.path(res_dir,"report_rmd_v2.Rmd"),
                  output_file = file.path(res_dir,"report_fileName_nonIFN.html"))
file_ori <- file.path(res_dir,"var_cand_IFN.xlsx")
rmarkdown::render(file.path(res_dir,"report_rmd_v2.Rmd"),
                  output_file = file.path(res_dir,"report_fileName_IFN.html"))
##
## END
