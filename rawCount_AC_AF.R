## Revised counting algorithm 
##
## Counting
count_dir <- file.path(base_dir,"rawCount")
dir.create(count_dir)
site_tab <- data.frame(variant=names(rd),stringsAsFactors = FALSE)
site_tab$chr <- gsub("(.*):(.*)_(.*)/(.*)","\\1",site_tab$variant)
site_tab$posi <- gsub("(.*):(.*)_(.*)/(.*)","\\2",site_tab$variant)
site_tab$ref <- gsub("(.*):(.*)_(.*)/(.*)","\\3",site_tab$variant)
site_tab$alt <- gsub("(.*):(.*)_(.*)/(.*)","\\4",site_tab$variant)
site_tab$idx <- paste0(site_tab$chr,":",site_tab$posi)
write.table(site_tab[,c(2,3,3)],file.path(count_dir,"varSite.list"),
            row.names = FALSE,sep = "\t",col.names = FALSE,quote = FALSE)
readCount <- function(bam_file,var_list=NULL,ref_seq=NULL,res_dir=count_dir){
  samID <- gsub("_recal.bam","",basename(bam_file))
  system(paste("bam-readcount",
               "-f",ref_seq,"-l",var_list,bam_file,">",
               file.path(res_dir,paste0(samID,"_rawCount.txt")),sep=" "))}
#
library(BiocParallel)
library(batchtools)
temp_file1 <- "/rugpfs/fs0/brc/scratch/jluo/bt-simple_bigmem.tmpl"
param <- BatchtoolsParam(workers=10000, cluster="slurm",template=temp_file1)
register(param)
#
# Calculate raw counts for each position
bam_file <- dir(file.path(base_dir,"BAM"),pattern="_recal.bam$",full.names = TRUE)
bplapply(bam_file,readCount,res_dir=count_dir,
         var_list=file.path(count_dir,"varSite.list"),
         ref_seq="path to BSgenome.Hsapiens.UCSC.hg19.fa")
#
# Generate DP and AD information
count_file <- dir(file.path(base_dir,"rawCount"),pattern="_rawCount.txt",full.names = TRUE)
countProc <- function(count_file,res_dir=NULL,site=site_tab){
  samID <- gsub("_rawCount.txt","",basename(count_file))
  count_tab <- read.delim(count_file,header = FALSE,stringsAsFactors = FALSE)
  count_dat <- data.frame(idx=paste0(count_tab$V1,":",count_tab$V2),
                          DP=as.integer(count_tab$V4),
                          stringsAsFactors = FALSE)
  count_dat$A <- as.integer(unlist(sapply(strsplit(count_tab$V6,split = ":"),"[",2)))
  count_dat$T <- as.integer(unlist(sapply(strsplit(count_tab$V7,split = ":"),"[",2)))
  count_dat$G <- as.integer(unlist(sapply(strsplit(count_tab$V8,split = ":"),"[",2)))
  count_dat$C <- as.integer(unlist(sapply(strsplit(count_tab$V9,split = ":"),"[",2)))
  count_mg <- merge(site,count_dat,by='idx',all.x=TRUE)
  count_mgA <- count_mg[count_mg$alt=="A",]
  count_mgA$AD <- count_mgA$A
  count_mgA$AF <- count_mgA$AD/count_mgA$DP
  count_mgT <- count_mg[count_mg$alt=="T",]
  count_mgT$AD <- count_mgT$T
  count_mgT$AF <- count_mgT$AD/count_mgT$DP
  count_mgG <- count_mg[count_mg$alt=="G",]
  count_mgG$AD <- count_mgG$A
  count_mgG$AF <- count_mgG$AD/count_mgG$DP
  count_mgC <- count_mg[count_mg$alt=="C",]
  count_mgC$AD <- count_mgC$C
  count_mgC$AF <- count_mgC$AD/count_mgC$DP
  count_mg <- rbind(count_mgA,count_mgT,count_mgG,count_mgC)
  count_mg <- dplyr::select(count_mg,variant,DP,AD,AF)
  names(count_mg) <- c("variant",paste0("DP_",samID),paste0("AD_",samID),paste0("AF_",samID))
  count_mg[is.na(count_mg$DP),-1] <- rep(0,3)
  write.table(count_mg,file.path(res_dir,paste0(samID,"_countProc.txt")),
              sep="\t",row.names = FALSE,quote = FALSE)}
# dir.create(file.path(base_dir,"rawCount"))
bplapply(count_file,countProc,res_dir=file.path(base_dir,"rawCount"),site=site_tab)
dep_file <- dir(file.path(base_dir,"rawCount"),pattern="_countProc.txt",full.names = TRUE)
count_matrix <- site_tab
for(k in 1:length(dep_file)){
  dep_mat <- read.delim(dep_file[k],stringsAsFactors = FALSE)
  count_matrix <- merge(count_matrix,dep_mat,by='variant')}
names(count_matrix) <- gsub("\\.","-",names(count_matrix))
write.table(count_matrix,file.path(base_dir,"rawCount","var_depComb.txt"),
            sep="\t",quote=FALSE,row.names = FALSE)
# save(list = c("rd","count_matrix"),file = file.path(base_dir,"rawCount_depComb_2019.RData"))
##
## VCFs from genereaed by Mutect2
base_dir <- "path to main directory"
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
vcf <- readVcf(file.path(base_dir,"Merged_anno.hg19_multianno.vcf"),"hg19")
rd <- rowRanges(vcf)
rd$seq <- getSeq(Hsapiens,seqnames(rd),start=start(rd)-5,end=end(rd)+5)
#
matDP <- count_matrix
#
matAC <- lapply(geno(vcf)$AD,`[[`,2)
rownames(matAC) <- seqnames(rd)
#
matAF<- lapply(geno(vcf)$AF,`[[`,1)
rownames(matAF) <- seqnames(rd)
#
matDP <- matDP[rownames(matDP) %in% rownames(matAF),]
matAC <- matAC[rownames(matAC) %in% rownames(matDP),]
matAF <- matAC[rownames(matAF) %in% rownames(matDP),]
save(matDP,matAC,matAF,file="path to variant RData")
##
## END
