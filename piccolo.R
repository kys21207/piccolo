#piccolo.R
#' @title PICCOLO: Colocalization of GWAS and GWAS/xQTL signals without summary stats
#'
#' @description PICCOLO estimates the probability that two genetic signals between GWAS and GWAS/xQTLs 
#'   are shared using colocalization approach. A proper colocalization analysis requires 
#'   full summary stats for both the GWAS and QTL studies. Unfortunately, studies frequently 
#'   do not share full summary stats which makes a traditional colocalization analysis impossible..
#'   In this case, PICCOLO enable to do a colocalization using PICS.
#'      
#' @param rs SNP rsid. Single string or vector.
#' @param chrom SNP chromosome. Single string or vector. 
#' @param pos SNP position. Single string or vector.
#' @param pval the association p-value for the SNP. Single string or vector.
#' @param ancestry Ancestry information: EUR/AMR/AFR/EAS/SAS. 
#' @param indication Associated phenotype or trait info. OPTIONAL!!
#' @param nCORE number of cores for multiple tasks (recommand nCORE=20)
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL) 
#'   
#' @return  
#'  a data frame containing the result of the
#'  PICS calculation.
#' 
#' @imports glue
#' @imports tidyverse
#' @imports gtx
#' @imports foreach
#' @imports doParallel
#' @imports magrittr
#' @author Kijoung Song & Xiao Zhang \email{kys21207@@gsk.com}
#' @export 

piccolo1 <- function(chrom,pos,rs,pval,ancestry,indication,hg38 = TRUE,dbc=getOption("gtx.dbConnection", NULL), nCORE = 1){
  
  # check database connection
  gtxdbcheck(dbc)
  
  # check input
  if((missing(rs) & missing(chrom) & missing(pos)) | (missing(rs) & missing(chrom)) |(missing(rs) & missing(pos))) { stop("piccolo | must specify either: rs or chr/pos", call. = FALSE) }
  if(!missing(rs) & !missing(chrom) & !missing(pos))  { gtx_fatal_stop("piccolo | Duplicate input: input either rs or chr/pos", call. = FALSE) }
  if(missing(pval)) { gtx_fatal_stop("piccolo | input is missing: pval", call. = FALSE) }
  if(!missing(rs)){
    if(missing(ancestry)) {
      gtx_fatal_stop("piccolo | required input: input ancestry (EUR/AFR/AMR/EAS/SAS)", call. = FALSE)
	}else{
  	if(!all(toupper(ancestry) %in% c("EUR","AFR","AMR","EAS","SAS"))) stop("piccolo | Check input: ancestry (EUR/AFR/AMR/EAS/SAS)", call. = FALSE)
	}
    if(missing(indication)) {indication <- c(rep(NA_character_, length(rs)))}
	tmp.list <- list(rs,pval,ancestry,indication)
	
    if(all(sapply(tmp.list,length)==length(tmp.list[[1]]))){
  	snpID <- rs
  	input <- tibble(snpID, pval, ancestry, indication)
  	
	}else{
      gtx_fatal_stop("piccolo | the length of each input is different.", call. = FALSE)
	}
	input <- input %>% distinct(snpID, indication, ancestry, .keep_all = TRUE)
	
  }else if(!missing(chrom) & !missing(pos)){
    if(missing(ancestry)) { gtx_fatal_stop("piccolo | required input: input ancestry (EUR/AFR/AMR/EAS/SAS)", call. = FALSE) }
    if(missing(indication)) { indication <- c(rep(NA_character_, length(pos)))}
	tmp.list <- list(chrom,pos,pval,ancestry,indication)
    
    if(all(sapply(tmp.list,length)==length(tmp.list[[1]]))){
  	snpID <- paste(chrom,pos,sep=":")
  	input <- tibble(snpID, pval, ancestry, indication)
	}else{
      gtx_fatal_stop("piccolo | the length of each input is different.", call. = FALSE)
	}
	input <- input %>% distinct(snpID, indication, ancestry, .keep_all = TRUE) 
  }
 
  input <- input %>% filter(!is.na(snpID) & !is.na(pval) & !is.na(ancestry)) %>%
                     mutate(ancestry=toupper(ancestry))
 
  input.check4rs <- input %>% filter(grepl("rs",snpID)) %>%
	                            mutate(snpID = gsub("[^[:alnum:][:blank:]+?&/\\-].*","",snpID))
  input.check4chrpos <- input  %>% filter(!grepl("rs", snpID) & grepl(":",snpID))
  excluded.snpID <- input %>% filter(!(snpID %in% input.check4rs$snpID) & !(snpID %in% input.check4chrpos$snpID))
 
  input <- rbind(input.check4rs,input.check4chrpos)
 
  # now seperate ancestry and run them individually
  # --------------------------------------------------------------------------------------------------
  final_res_temp <- list()
  final_n <- 1
  for(i in unique(input$ancestry)){
    input_single_ancestry <- input %>% filter(ancestry == i)
	temp_res <- piccolo_single_ancestry(input_single_ancestry, hg38 = hg38, nCORE = nCORE)
 
    if(!is.null(temp_res)) {
      final_res_temp[[final_n]] <- temp_res
   	final_n <- final_n + 1
	}
   }
 
  final_res <- int_fastDoCall("rbind", final_res_temp)
      
  # I add
  return(final_res)
} 
 

#' Nested in piccolo
#'
#'     Estimate colocalization of GWAS and GWAS/eQTL/pQTL signals based on single ancestry
#'     Nested in piccolo
#' 
#' @param input tibble(snpID, pval, ancestry, indication), snpID=paste(chrom,pos,sep=":")
#' @param nCORE number of cores
#' @return 
#'   H4 based on each ancestry 
#' @author  Kijoung Song & Xiao Zhang \email{kys21207@@gsk.com} 
#' 
#' @export


piccolo_single_ancestry <- function(input,hg38 = TRUE, dbc=getOption("gtx.dbConnection", NULL), nCORE = 1){
 
 
# check in single ancestry
# keep it as a seperate variable
if(length(unique(input$ancestry)) != 1){stop("piccolo_single_ancestry | this function require single ancestry")}
 
 
#-----------------------------------------------------------------
# calculate pics by using pics_calc1
# In the last step, we already checked input data has unique ancestry
 
gwas.pics <- tryCatch(pics_calc1(input, hg38 = hg38, nCORE_pics = nCORE),error=function(e) NULL)
 
if(is.null(gwas.pics)) {
 gtx_warn("All SNP IDs could not be found in the current LD reference.\n
            Check your input data carefully. \n
            For example, either hg37 or hg38, ancestry(=EUR/AMR/AFR/EAS/SAS)")
return(NULL)
}
 
  
gwas.pics <- gwas.pics %>% mutate(ID=paste(chrom_idx,pos_idx,pval_idx,indication,rsid_idx,snpID_idx,ancestry,sep = ";"))
#-----------------------------------------------------------------
# Find the nearest gene
# select SNP that are in LD with index SNP, group by ID, then get min and max pos  	
gwas.pics.nest <- gwas.pics %>% select(ID,pos,pics) %>% group_by(ID) %>% nest()
idx.snps <- gwas.pics %>% group_by(chrom_idx, pos_idx, pval_idx, indication, rsid_idx, snpID_idx, ancestry) %>%
                          summarise(start = min(pos), end = max(pos))
 
#-------------------------
# find ensemblid that are near index SNP (+-500 bandwidth from min and max pos) from genes table
 
# genes table contains position info and ensemble ID
# use "protein_coding" as genetype
# idx.nearst.genes contains ensembleid and ID (a combination of chrom_idx, pos_idx, rsid_idx, snpID_idx)
 
 
idx.nearst.genes <- idx.snps %>%
  group_by(chrom_idx, pos_idx, pval_idx, indication, rsid_idx, snpID_idx,ancestry) %>%
do(sqlWrapper(dbc, paste0("SELECT ensemblid FROM gene_gwas_hg38_use.t_genes WHERE genetype = 'protein_coding' AND chrom = '",
                            .$chrom_idx, "' AND (", "(pos_start >= ", .$start - 5e+05,
                            " AND pos_start <= ", .$end + 5e+05, ")", " OR (pos_end >= ",
                            .$start - 5e+05, " AND pos_end <= ", .$end + 5e+05, ")",
                            " OR (pos_start >= ", .$start - 5e+05, " AND pos_end <= ",
                            .$end + 5e+05, ")", ")"), uniq = F, zrok = TRUE)) %>%
  mutate(ID=paste(chrom_idx,pos_idx,pval_idx,indication, rsid_idx, snpID_idx, ancestry,sep = ";"))
 
#  k_pics is ancestry-specified pics table
#  pics table is xQTL data constaining pics, index SNP (useless), pos etc. info for each entity (i.e. ensemble ID)
#  look for ensemble ID from that is in k_pics
tmp <- NULL
pics.qtls<-NULL 
 
if (length(unique(idx.nearst.genes$ensemblid)) > 5000) {
ens2List <- split(unique(idx.nearst.genes$ensemblid),
                    cut(1:length(unique(idx.nearst.genes$ensemblid)), ceiling(length(unique(idx.nearst.genes$ensemblid))/5000),F))
} else {
ens2List <- split(unique(idx.nearst.genes$ensemblid),1)
}
 

if(hg38){
 
# PRL
registerDoParallel(cores=nCORE)
 
pics.qtls <- foreach(j=1:length(ens2List), .combine = rbind) %dopar% {
tmp <- sqlWrapper(dbc, paste0("SELECT * FROM gene_gwas_hg38_use.pics_qtls_global WHERE ensemblid IN ('",
                                paste(ens2List[[j]], collapse = "', '"), "') and ethnicity='", unique(input$ancestry), "'"),
                    uniq = F, zrok = TRUE) %>% mutate(ancestry = unique(input$ancestry))
#tmp
}
stopImplicitCluster()

pics.qtls <- pics.qtls %>% mutate(ID=paste(ancestry,analysis,ensemblid,rsid_index,pval_index,pubmedid,hgncid,entity,xqtl_type, sep = ";"))

  #head(pics.qtls)
# for each tissue/entity we keep only the highest pics and then stores in pos.idx
pos.idx <- pics.qtls %>%arrange(ID,-pics) %>% group_by(ID)%>% filter(row_number()==1) %>% select(ID,pos)
#pos.idx
                         
 
# n1 is the number of ID in xQTL
# idx.nearst.genes has the nearest genes (500 by pos) of each entity/tissue from genes table
# gwas.pics.nest has ID and the calculated pics from index.data (a GWAS data)
 
idx.nearst.genes.id <- unique(idx.nearst.genes$ID)
# PRL
registerDoParallel(cores=nCORE)
 
res <- foreach(j=1:length(idx.nearst.genes.id), .combine = rbind) %dopar% {
     dta.gwas <- subset(gwas.pics.nest, ID %in% idx.nearst.genes.id[j])
     tmp <- subset(idx.nearst.genes, ID %in% idx.nearst.genes.id[j])
	
	# limit entity in pics.qtl (= pics table = xQTL) to that is nearest gene to index SNP (GWAS data)  
  dta.qtls <- pics.qtls %>% filter(ensemblid %in% tmp$ensemblid) %>% select(ID,pos,pics) %>%
            group_by(ID) %>% nest() %>% rename(ID1=ID,data1=data)
	
	# crossing is looking for all possible combinations of GWAS and xQTL
	# GWAS data is gwas.pics.nest, QTL data is idx.nearst.genes
tmp00 <- crossing(dta.gwas, dta.qtls) %>%
	mutate(H4 = map2_dbl(data, data1, int_pico_lite)) %>% select(ID,ID1,H4) %>% distinct()
#tmp00
}
stopImplicitCluster()                   
 
  
}else{
  
  type = case_when(
    unique(input$ancestry) == 'EUR' ~ "pics_qtls",
    unique(input$ancestry) == 'AFR' ~ "pics_afr",
    unique(input$ancestry) == 'AMR' ~ "pics_amr",
    unique(input$ancestry) == 'SAS' ~ "pics_sas",
    unique(input$ancestry) == 'EAS' ~ "pics_eas"
  )
 # PRL
registerDoParallel(cores=nCORE)
 
pics.qtls <- foreach(j=1:length(ens2List), .combine = rbind) %dopar% {
tmp <- sqlWrapper(dbc, paste0("SELECT * FROM gene_gwas_use.",type," WHERE entity IN ('",
                                paste(ens2List[[j]], collapse = "', '"), "')"),
                    uniq = F, zrok = TRUE) %>% mutate(ancestry = unique(input$ancestry),xqtl_type=NA)
#tmp
}
stopImplicitCluster()
                      
pics.qtls <- pics.qtls %>% mutate(ID=paste(ancestry,tissue,entity,rsid_idx,pval_idx,pmid_idx,hgncid,entity,xqtl_type, sep = ";"))
#head(pics.qtls)
# for each tissue/entity we keep only the highest pics and then stores in pos.idx
pos.idx <- pics.qtls %>% arrange(ID,-pics) %>% group_by(ID)%>% filter(row_number()==1) %>% select(ID,pos)
#pos.idx
                         
 
# n1 is the number of ID in xQTL
# idx.nearst.genes has the nearest genes (500 by pos) of each entity/tissue from genes table
# gwas.pics.nest has ID and the calculated pics from index.data (a GWAS data)

idx.nearst.genes.id <- unique(idx.nearst.genes$ID)
# PRL
registerDoParallel(cores=nCORE)
 
res <- foreach(j=1:length(idx.nearst.genes.id), .combine = rbind) %dopar% {
     dta.gwas <- subset(gwas.pics.nest, ID %in% idx.nearst.genes.id[j])
     tmp <- subset(idx.nearst.genes, ID %in% idx.nearst.genes.id[j])
	
	# limit entity in pics.qtl (= pics table = xQTL) to that is nearest gene to index SNP (GWAS data)  
  dta.qtls <- pics.qtls %>% filter(entity %in% tmp$ensemblid) %>% select(ID,pos,pics) %>%
            group_by(ID) %>% nest() %>% rename(ID1=ID,data1=data)
	
	# crossing is looking for all possible combinations of GWAS and xQTL
	# GWAS data is gwas.pics.nest, QTL data is idx.nearst.genes
tmp00 <- crossing(dta.gwas, dta.qtls) %>%
         mutate(H4 = map2_dbl(data, data1, int_pico_lite)) %>% select(ID,ID1,H4) %>% distinct()
##tmp00
}
stopImplicitCluster()                   

}
                                   
#----------------------------
# split ID, this is in gwas data; split ID1, this is in xQTL data   
res <- res %>% inner_join(pos.idx,by=c("ID1"="ID")) %>% rename(qtl_pos=pos) %>%
  separate(ID,c("gwas_chrom","gwas_pos","gwas_pval","gwas_analysis","gwas_rsid","snpID","ancestry"),';') %>%
  separate(ID1,c("ancestry1","tissue","qtl_ensemblid","qtl_rsid","qtl_pval","qtl_pubmedid","qtl_hgncid","entity","xqtl_type"),';') %>%
  select(snpID,gwas_analysis,gwas_chrom,gwas_rsid,gwas_pos,gwas_pval,tissue,qtl_rsid,qtl_pos,qtl_pval,
          qtl_hgncid,qtl_ensemblid,qtl_pubmedid,ancestry,entity,xqtl_type,H4) %>%
  rename(gwas_input=snpID) %>%
  mutate(gwas_pos=as.numeric(gwas_pos),gwas_pval=as.numeric(gwas_pval),qtl_pos=as.numeric(qtl_pos),qtl_pval=as.numeric(qtl_pval))
  
# display SNPs that exist in xQTL but not in GWAS
check.missing.gwas.snp <- input %>% filter(!(snpID %in% unique(res$gwas_input))) %>%
                           select(snpID,pval,ancestry,indication) %>%
                           rename(gwas_input=snpID,gwas_pval=pval,gwas_analysis=indication)
 
if(nrow(check.missing.gwas.snp) >= 1) res <- int_sbind(check.missing.gwas.snp,res)
 gtx_info('PICCOLO analysis is complete.')
return(res)

}
                        
                        
#==========================================================================================================================
#' @title PICS calculation using the linkage information () 
#' 
#' @description  The PICS algorithm calculates the most likely causal SNPs given the observed association signal
#'   at a locus. For an associated locus, enter the most highly-associated SNP 
#'   (referred to as the index SNP) and the strength of association. Using the linkage information
#'   in 10,000 UKB samples for EUR and in 1000 genome reference P3 for Non-EUR, the algorithm identifies 
#'   the SNPs that are most likely to be the causal variants responsible for the association (PICS_Probability).
#'
#'   See \strong{Genetic and Epigenetic Fine-Mapping of Causal Variants in Autoimmune Disease} 
#'   by Kyle Kai-How Farh,et al. Nature 518, 337–343 (19 February 2015) at
#'   \url{http://www.nature.com/nature/journal/v518/n7539/full/nature13835.html#close}
#' @param index.data TBD
#' @param dbc Database connection. Default: getOption("gtx.dbConnection", NULL) 
#' @param nCORE_pics Define the number of cores for multiple tasks
#'   
#' @return  
#'  a data frame containing the result of the
#'  PICS calculation.
#' 
#' @author Kijoung Song & Xiao Zhang \email{kys21207@@gsk.com}
#' 
#' @export
                        
pics_calc1 <- function(index.data,hg38 = TRUE, dbc=getOption("gtx.dbConnection", NULL), nCORE_pics = 1){
 
# case 1, when index data only containts "rs"
# -------------------------------------------------------------------------
rs.snpid <- subset(index.data, grepl("rs", index.data$snpID))
rs.snpid.tmp <- gsub("rs", "", rs.snpid$snpID)

dta.ext <- NULL  
 if (length(rs.snpid.tmp) > 0) {
  if (length(unique(rs.snpid.tmp)) > 5000) {
    rsList <- split(unique(rs.snpid.tmp), cut(1:length(unique(rs.snpid.tmp)), 
                                              ceiling(length(unique(rs.snpid.tmp))/5000), F))
  } else {
    rsList <- split(unique(rs.snpid.tmp), 1)
  }     

if(hg38){
registerDoParallel(cores=nCORE_pics)
dta.ext <- foreach(i=1:length(rsList), .combine = rbind) %dopar% {
    x <- tryCatch(sqlWrapper(dbc, paste0("SELECT chrom,pos,rs FROM gene_gwas_hg38_use.gwas_dbsnp WHERE rs IN (", 
                                         paste(rsList[[i]], collapse = ","), ")"), uniq = F, 
                             zrok = FALSE), error = function(e) NULL)
    if(is.null(x)){return(NULL)}

    if (!is.null(x)) {
      x <- x %>% filter(!is.na(rs)) %>% distinct(chrom, pos, rs) %>%
                 mutate(rsid1 = paste0("rs",rs), snpID = rsid1) %>%
                 select(-rs) 
    } 
}
stopImplicitCluster()  
}else{
registerDoParallel(cores=nCORE_pics)
dta.ext <- foreach(i=1:length(rsList), .combine = rbind) %dopar% {
    x <- tryCatch(sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM gene_gwas_use.sites WHERE rsid IN (", 
                                         paste(rsList[[i]], collapse = ","), ")"), uniq = F, 
                             zrok = FALSE), error = function(e) NULL)
    if(is.null(x)){return(NULL)}

    if (!is.null(x)) {
      x <- x %>% filter(!is.na(rsid)) %>% distinct(chrom, pos, rsid) %>%
                 mutate(rsid1 = paste0("rs",rsid), snpID = rsid1) %>%
                 select(-rsid) 
    } 
}
stopImplicitCluster()
}
}
# -------------------------------------------------------------------------

# rs.snpid.tmp and cp.snpid contains snpID, the numbers without "rs" for index SNPs
# case 2),  
cp.snpid <- subset(index.data, !grepl("rs", index.data$snpID))
dta.cp <- NULL                                 
#------------------------- 
# get snpID, chrom, and pos for index SNP by using site table
# finally save to dta.ext, then merge back to index.data
if(nrow(cp.snpid)>0){
# case 2), looks for rs for chrom/pos  
cp.snpid$chrom <- unlist(lapply(strsplit(as.character(cp.snpid$snpID), ":"), function(x) x[[1]]))
cp.snpid$pos <- as.numeric(unlist(lapply(strsplit(as.character(cp.snpid$snpID), ":"), function(x) x[[2]])))
 
tmp00 <- list() 
n <- 1  

for (i in unique(cp.snpid$chrom)) {
cp.snpid.sub <- subset(cp.snpid, chrom == i)
	
# rsList[(1-5000), (5001,10000),...], rsList[[1]] is the first 5000 pos within each chrom
if(length(unique(cp.snpid.sub$pos)) > 5000) {
  rsList <- split(unique(cp.snpid.sub$pos), cut(1:length(unique(cp.snpid.sub$pos)), ceiling(length(unique(cp.snpid.sub$pos))/5000), F))
}else{
  rsList <- split(unique(cp.snpid.sub$pos), 1)
}

if(hg38){
registerDoParallel(cores=nCORE_pics)
tmp01 <- foreach(j=1:length(rsList), .combine = rbind) %dopar% {

      x = sqlWrapper(dbc, paste0("SELECT chrom,pos,rs FROM gene_gwas_hg38_use.gwas_dbsnp WHERE chrom = '",
                                  	i, "' AND pos IN (", paste(rsList[[j]], collapse = ","),
                                      ")"), uniq = F, zrok = TRUE)
	
      if(is.null(x)){return(NULL)}

# in case there is no chrom:pos for a rs snpID, add chrom:pos
# rsid1 has "rs" and the number	
if (!is.null(x)) {
	 x <- x %>% filter(!is.na(rs)) %>% distinct(chrom, pos, rs) %>%
              mutate(snpID = paste(chrom, pos, sep = ":"),
         	           rsid1 = paste0("rs", rs)) %>%
              select(-rs)
#	 x
}
}
stopImplicitCluster()
  
}else{
# PRL
registerDoParallel(cores=nCORE_pics)
tmp01 <- foreach(j=1:length(rsList), .combine = rbind) %dopar% {

      x = sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM gene_gwas_use.sites WHERE chrom = '",
                                  	i, "' AND pos IN (", paste(rsList[[j]], collapse = ","),
                                      ")"), uniq = F, zrok = TRUE)
	
      if(is.null(x)){return(NULL)}

# in case there is no chrom:pos for a rs snpID, add chrom:pos
# rsid1 has "rs" and the number	
if (!is.null(x)) {
	 x <- x %>% filter(!is.na(rsid)) %>% distinct(chrom, pos, rsid) %>%
              mutate(snpID = paste(chrom, pos, sep = ":"),
         	           rsid1 = paste0("rs", rsid)) %>%
              select(-rsid)
#	 x
}
}
stopImplicitCluster()
}
tmp00[[n]] <- tmp01
n <- n + 1
}
# res contains all 3 info no matter case 1 or case 2                
dta.cp <- int_fastDoCall("rbind", tmp00)                           	

# because we deleted case 1, here add this back
# dta.ext <- NULL                    	
if (is.null(dta.cp)) gtx_warn("pics_calc: all chr:pos' are missing, check your chr:pos' carefully!")
}
  
if (is.null(dta.ext) & !is.null(dta.cp)) {
   dta.ext <- dta.cp
} else if(!is.null(dta.ext) & !is.null(dta.cp)){
# finally saved all 3 info into dta.ext 
   dta.ext <- rbind(dta.ext, dta.cp)
}
#dta.ext

# add back to index.data. Now we have rsid1                 
if (!is.null(dta.ext)) {
  index.data = merge(index.data, dta.ext, by = "snpID")
} else {
  gtx_fatal_stop("pics_calc: all rs's are missing, carefully check the input data, specially rs.", call. = FALSE)
}
#index.data

# exclude chrom "X" and "XY"
index.data <- index.data %>% filter(chrom %in% as.character(seq(1,22)))

#==========================================================================
# find chrom/pos that are in LD with index SNP
# LD table contains all r between pos 1 and pos 2
# duplicate chrom and pos for index SNP: chrom=chrom1=chrom2, pos=pos1=pos2
tmp.ld <- index.data %>% mutate(chrom2 = chrom, pos2 = pos, r=1, r2=1)%>%
          select(snpID,rsid1,pval,chrom,pos,chrom2,pos2,ancestry,indication,r,r2) %>%
          rename(chrom1=chrom,pos1=pos)

	
tmp00 <- list()
n <- 1 

k<- unique(index.data$ancestry)
# check again
if(length(k) != 1){gtx_fatal_stop("pics_calc | index.data should have unique ancestry")}

# pos1List stores chrom/pos that are from the index.data but in each ancestry
for (i in unique(index.data$chrom)) {
sub.dta <- index.data %>% filter(ancestry %in% k & chrom %in% i)
if (length(unique(sub.dta$pos)) > 5000) {
	pos1List <- split(unique(sub.dta$pos), cut(1:length(unique(sub.dta$pos)), ceiling(length(unique(sub.dta$pos))/5000), F))
} else {
	pos1List <- split(unique(sub.dta$pos), 1)
}

if(hg38){
# PRL
# sub.ld contains all pos2 for pos1 (i.e. index SNP)
registerDoParallel(cores=nCORE_pics)
tmp01 <- foreach(j=1:length(pos1List), .combine = rbind) %dopar% {
	sub.ld = sqlWrapper(dbc, paste0("SELECT * FROM gene_gwas_hg38_use.ld_global WHERE pos1 IN (", paste(pos1List[[j]], collapse = ","),
                                	") and chrom1='", i, "' and ethnicity='", k, "'"), uniq = F, zrok = TRUE)
    if(is.null(sub.ld)) {return(NULL)}
    sub.ld[!duplicated(sub.ld[, c("pos2", "pos1")]),]
#	sub.ld
}
stopImplicitCluster()

}else{

 type = case_when(
    k == 'EUR' ~ "ld",
    k == 'AFR' ~ "ld_afr",
    k == 'AMR' ~ "ld_amr",
    k == 'SAS' ~ "ld_sas",
    k == 'EAS' ~ "ld_eas"
  )
registerDoParallel(cores=nCORE_pics)
tmp01 <- foreach(j=1:length(pos1List), .combine = rbind) %dopar% {
	sub.ld = sqlWrapper(dbc, paste0("SELECT * FROM gene_gwas_use.",type," WHERE pos1 IN (", paste(pos1List[[j]], collapse = ","),
                                	") and chrom1='", i, "'"), uniq = F, zrok = TRUE)
    if(is.null(sub.ld)) {return(NULL)}
    sub.ld[!duplicated(sub.ld[, c("pos2", "pos1")]),]
#	sub.ld
}
stopImplicitCluster()

}

tmp00[[n]] <- tmp01
n <- n+1
}
#  }
all.ld <- int_fastDoCall("rbind", tmp00)

#----------------------------------- 
all.ld <- all.ld %>% mutate(r2 = r^2) %>% arrange(chrom1, pos1,-r2)

# after merging with index.data, all.ld has pos1 with its all pos2 (including itself), and r2, snpID	
all.ld <- all.ld %>% inner_join(index.data, by = c("chrom1"="chrom", "pos1"="pos")) 
all.ld <- int_sbind(tmp.ld, all.ld)

# those with r2>0.5 are in LD with index SNP              	
all.ld <- all.ld %>% filter(r2 > 0.5) %>%
                     mutate(pval = ifelse(pval == 0, 1e-250, pval)) %>%
                     mutate(SD = ifelse(r > 0.99999,0,sqrt(1 - abs(r)^6.4) * sqrt(-log10(pval))/2)) %>%
                     mutate(Mean = r2 * (-log10(pval))) %>%
                     mutate(prob = ifelse(SD == 0, 0.8, 1 - pnorm(-log10(pval), Mean, SD))) 
prob.sum <- all.ld %>% group_by(chrom1, pos1, pval) %>% summarise(prob_sum = sum(prob))
all.ld <- all.ld %>% inner_join(prob.sum, by = c("chrom1", "pos1", "pval")) %>%
              # pics = prob/prob_sum, prob_sum is the sum of all SNP in LD
              # the sum of pics across SNPs sharing the same index SNP is 1
                     mutate(pics = prob/prob_sum)
 
# add rs back to all.ld that contains chrom and pos only
# complete info of SNP in LD with index SNP
tmp00 <- list()
n <- 1
 
for (i in unique(all.ld$chrom2)) {
sub.dta <- subset(all.ld, chrom2 %in% i)
if (length(unique(sub.dta$pos2)) > 5000) {
	pos2List <- split(unique(sub.dta$pos2), cut(1:length(unique(sub.dta$pos2)), ceiling(length(unique(sub.dta$pos2))/5000), F))
} else {
	pos2List <- split(unique(sub.dta$pos2), 1)
}
tmp01 <- list()    	
for (j in 1:length(pos2List)) {
if(hg38){
  tmp <- sqlWrapper(dbc, paste0("SELECT chrom,pos,rs FROM gene_gwas_hg38_use.gwas_dbsnp WHERE pos IN (",
                                  paste(pos2List[[j]], collapse = ","), ") and chrom='",
                              	i, "'"), uniq = F, zrok = FALSE) %>% rename(rsid=rs)
}else{
  tmp <- sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM gene_gwas_use.sites WHERE pos IN (",
                                  paste(pos2List[[j]], collapse = ","), ") and chrom='",
                              	i, "'"), uniq = F, zrok = FALSE)
}
#  	tmp <- sqlWrapper(dbc, paste0("SELECT chrom,pos,rsid FROM sites_onestop_annot WHERE pos IN (",
#                                    paste(pos2List[[j]], collapse = ","), ") and chrom='",
#                                	i, "'"), uniq = F, zrok = FALSE)
	tmp1 <- tmp %>% filter(!is.na(rsid)) %>% distinct(chrom, pos, rsid) %>%
	                mutate(rsid=paste0("rs",rsid))
	tmp2 <- tmp %>% filter(is.na(rsid)) %>%
	                filter(!(paste(chrom, pos, sep = "_") %in% paste(tmp1$chrom, tmp1$pos, sep = "_"))) %>%
                  mutate(rsid=paste0("rs",rsid))
	tmp01[[j]] <- rbind(tmp1, tmp2)
}
tmp00[[n]] <- int_fastDoCall("rbind", tmp01)
n <- n+1
}
snp.anno <- int_fastDoCall("rbind", tmp00)
# now snp.anno comtains all information for chrom2, pos2, rsid, and snp2 for a later binding with xQTL         	
snp.anno <- snp.anno %>% mutate(snp2 = rsid) %>%
                         rename(chrom2=chrom,pos2=pos)
 
#-------------------------
# pics.result renamed XX2 to XX, and XX to XX_idx
# each line is different: chrom, pos, rsid (rsXXX), pics, ancestry, indication
# share one: snpID_idx, pval_idx, chrom_idx, pos_idx
pics.result <- all.ld %>% inner_join(snp.anno, by = c("chrom2", "pos2")) %>% arrange(snpID, -r2) %>%
                          select(chrom2, pos2, snp2, pics, ancestry, indication, 
                                 snpID, rsid1, pval, chrom1, pos1) %>%
                          rename(chrom=chrom2,pos=pos2,rsid=snp2,snpID_idx=snpID,
                                rsid_idx=rsid1,pval_idx=pval,chrom_idx=chrom1,pos_idx=pos1)
#  pics.result
 
return(pics.result)
}
 
    
#=============================================================================================================================  
#' int_coloc_pics_lite 
#'   Test for colocalization of two PICS sets
#' 
#' @return only H3 & H4 posteriors
#' @param data1,data2  PICS sets from read.pics or download.pics
#' @param pics1,pics2  column name to pull PICS prob from. Default = "PICS_probability"
#' @param rsid1,rsid2  column name to pull rsid from.      Default = "Linked_SNP"
#' @param rounded   Decimal points to round posteriors to
#' @param priorc1   Prior probability for colocalization with siganl for data1  Default = 1e-4
#' @param priorc2   Prior probability for colocalization with siganl for data2 Default = 1e-4
#' @param priorc12  Prior probability for colocalization of both signals.   Default = 1e-5
#' 
#' @author Karsten Sieber \email{karsten.b.sieber@@gsk.com}  


int_coloc_pics_lite1 <- function(data1,
                                data2,
                                pics1    = "PICS_probability", # column header for poster probabilities in data1
                                pics2    = "PICS_probability", # column header for poster probabilities in data2
                                rsid1    = "Linked_SNP",       # column header for snps in LD in data1
                                rsid2    = "Linked_SNP",       # column header for snps in LD in data2
                                rounded  = 6,
                                priorc1  = 1e-4, 
                                priorc2  = 1e-4, 
                                priorc12 = 1e-5
) {
  stopifnot(exists("data1") & exists("data2"))
  if(is.logical(data1)){
    if(is.na(data1)){
      gtx_warn("int_coloc_pics_lite: data1 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  if(is.logical(data2)){
    if(is.na(data2)){
      gtx_warn("int_coloc_pics_lite: data2 is NA, skipping coloc.")
      return(list(results = NA, nvariants = NA))
    }
  }
  pics <- int_harmonize_pics(data1, 
                             data2, 
                             opts <- data.frame(rsid1 = rsid1, rsid2 = rsid2, pics1 = pics1, pics2 = pics2, stringsAsFactors = FALSE))
  
  nv <- dim(pics)[1]
  res <- data.frame(prior = norm1(c(priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf    = c((sum(pics[[1]])*sum(pics[[2]]) - sum(pics[[1]]*pics[[2]]))/(nv*(nv - 1)), 
                              sum(pics[[1]]*pics[[2]])/nv))
  res$bf <- res$bf/res$bf[1]
  res$posterior <- norm1(res$prior*res$bf)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(res$posterior[2])
}

#' int_pico_lite 
#'     Simplify int_coloc_pics_lite function so that we can use map2
int_pico_lite <- function(x,y) int_coloc_pics_lite1(x,y,pics1="pics",pics2="pics",rsid1="pos",rsid2="pos")
                     
  
#' @title int_harmonize_pics
#' @description  int_harmonize_pics
#' @param data1 TBD
#' @param data2 TBD
#' @param opts TBD
int_harmonize_pics <- function(data1,
                               data2, 
                               opts = data.frame(pics1 = "PICS_probability",
                                                 pics2 = "PICS_probability",
                                                 rsid1 = "Linked_SNP",
                                                 rsid2 = "Linked_SNP",
                                                 stringsAsFactors = FALSE)
){
  ids <- unique(c(data1[[opts$rsid1]], data2[[opts$rsid2]]))
  tmp <- as.data.frame(matrix(data = NA, nrow = length(ids), ncol = 2))
  pp1 <- if (opts$pics1==opts$pics2) paste(opts$pics1, ".1", sep = "") else opts$pics1
  pp2 <- if (opts$pics1==opts$pics2) paste(opts$pics2, ".2", sep = "") else opts$pics2
  colnames(tmp) <- c(pp1, pp2)
  for(n in 1:length(ids)){
    tmp[[pp1]][n] <- if(!is.na(match(ids[n], data1[[opts$rsid1]]))) data1[which(data1[[opts$rsid1]]==ids[n]),][[opts$pics1]][1] else 0
    tmp[[pp2]][n] <- if(!is.na(match(ids[n], data2[[opts$rsid2]]))) data2[which(data2[[opts$rsid2]]==ids[n]),][[opts$pics2]][1] else 0 
  }
  res <- as.data.frame(cbind(
    norm1(tmp[[pp1]], log = FALSE),
    norm1(tmp[[pp2]], log = FALSE)
  ))
  colnames(res) <- c(pp1, pp2)
  rownames(res) <- ids
  return(res)
}



#' Combine Shape Objects
#'
#'   Combine shape objects into one shape object. It works analogous to rbind.
#'
#' @author  Martijn Tennekes
#'

#' @title int_sbind
#' @description int_sbind
#' @param x TBD
#' @param y TBD
#' @param fill TBD
int_sbind = function(x, y, fill=NA) {
  int_sbind.fill = function(d, cols){ 
    for(c in cols)
      d[[c]] = fill
    d
  }
  
  x = int_sbind.fill(x, setdiff(names(y),names(x)))
  y = int_sbind.fill(y, setdiff(names(x),names(y)))
  
  rbind(x, y)
}

#' An Alternative To The Internal Do.Call
#'  
#'   The do.call can be somewhat slow, especially when working with large objects. 
#'   This function is based upon the suggestions from Hadley Wickham on the R mailing list, 
#'   see here. Also thanks to Tommy at StackOverflow for suggesting how to handle double 
#' . and triple colon operators, ::, further enhancing the function.
#'  
#' @param what   either a function or a non-empty character string naming the function to be called.
#' @param args   a list of arguments to the function call. The names attribute of args gives the argument names.
#' @param quote  a logical value indicating whether to quote the arguments.
#' @param envir  an environment within which to evaluate the call. This will be most useful if what is a character string and the arguments are symbols or quoted expressions.
#'
#' @author  Max Gorden 
#' 
#

int_fastDoCall <- function(what, args, quote = FALSE, envir = parent.frame()){
  if (quote)
    args <- lapply(args, enquote)
  
  if (is.null(names(args))){
    argn <- args
    args <- list()
  }else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }
  
  if (class(what) == "character"){
    if(is.character(what)){
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if(length(fn)==1) {
        get(fn[[1]], envir=envir, mode="function")
      } else {
        get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
      }
    }
    call <- as.call(c(list(what), argn))
  }else if (class(what) == "function"){
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  }else if (class(what) == "name"){
    call <- as.call(c(list(what, argn)))
  }
  
  eval(call,
       envir = args,
       enclos = envir)
}
  
#' norm1
#' @author: unknown
  norm1 <- function(x, log = FALSE) {
	if (all(is.na(x))) return(x)
	if (log) {
  	x <- x - max(x, na.rm = TRUE)
  	x <- exp(x)
	} else {
  	
  	stopifnot(all(x >= 0, na.rm = TRUE))
  	x <- x / max(x, na.rm = TRUE)
	}
	return(x / sum(x, na.rm = TRUE))}
	
