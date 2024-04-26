# gtx
Genetics ToolboX R package
# PICCOLO: Colocalization of GWAS and xQTL signals without summary stats

## Guide for PICCOLO on CDSW
PICCOLO is a Rscript that can be run on CDSW. PICCOLO will take simple rsIDs or chr/pos and pvalues and determine the genetic signal(s) colocalization with any xQTLs.  
The basic workflow of the script is as follows.

For each genetic loci (input rsID and pvalue) determine the PICS credible set
Find all genes +/- 500 MB from the credible set
For each gene, find all xQTLs in each tissue
For each tissue with an xQTL for each gene, colocalize the input genetic signal and the xQTL

## Run PICCOLO
To run PICCOLO, use m37A and m37B modules or following command and input options.
```R
piccolo(rs="rs1234",pval=1.0e-9,indication="test",ancestry="EUR",hg38=TRUE, nCORE=7)  
```
1. indication, hg38 and nCORE are optional. 
2. ancestry requires one of 5 populations
​EUR, European
AFR, African
AMR, Ad Mixed American
EAS, East Asian
SAS, South Asian​ 
3. hg38=TRUE is the default. When working with hg19 data, then hg38=FALSE is required. 
Setting up nCORE requires first to find out how many cores we want to use from the ones we have available. It is recommended to leave one free core for other tasks.
```R
parallel::detectCores()

[1] 8
```
​​nCORE=7

OR
```R
piccolo(chr="1", pos=1234,pval=1.0e-9,ancestry="AFR",nCORE=10)
```
### input allows single value/string or vector 
​Run time is highly variable based on the number of genes and xQTLs across tissues in each locus.

## PICCOLO output
For every input locus, for every gene and tissue combination, a new row will be output for the colocalization analysis. The six output columns are:

* gwas_analysis: input indication of gwas
* gwas_chrom: input chrom of gwas
* gwas_rsid: index SNP ID for gwas  
* gwas_pos: input position of gwas
* gwas_pval: input p-value of gwas
* qtl_analysis: Tissue abbreviation for the gene xQTL tested. Full tissue key is here: QTLs key and info.xlsx​
* qtl_rsid: Index snp for the xQTL in the given tissue
* qtl_pos: position for the qtl_rsid & tissue
* qtl_pval: p-value for the the qtl_rsid & tissue
* qtl_hgncid: The gene symbol for the xQTL based on the qtl_irsid& tissue
* qtl_ensemblid: ensembl ID for the qtl_rsid & tissue​
* qtl_pubmedid: pubmed ID for the qtl data
* ancestry: input ancestry of interest
* H4: The posterior probability that the GWAS and eQTL signals DO colocalize

## Interpreting the output
1. Range for H4= 0-1
2. 0 = 0% probability
3. 1 = 100% probability
4. H4 cutoff: Colocalizations using full summary stats have found that H4>0.8 has been reasonable cutoffs. 
	We are currently performing analyses to determine if this is a proper cutoff when using PICCOLO.
Tip - The CSV output can be opened in EXCEL which allows for easy filtering of results.

