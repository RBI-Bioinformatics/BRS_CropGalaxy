#!/usr/local/bin/Rscript
library(getopt)

BY_CHROMOSOME_PLOTS = FALSE

if(!interactive()){
  opt = getopt(spec = matrix(c(
    "trait_name",  "t", 2, "character",
    "gwas_result", "f", 1, "character",
    "genomewide", "g", 2, "numeric",
    "suggestive", "u", 2, "numeric"
    , "bychrom", "c", 2, "logical"
  ), byrow = TRUE, ncol=4))
  
  # Parse options
  trait_name = opt$trait_name
  if(is.null(trait_name)){
    trait_name = paste0("Trait-", date())
  }
  gwas_result_file = opt$gwas_result
  if(is.null(gwas_result_file)){
    stop("Please specify GWAS result file as returned by GEMMA (columns chr,rs,ps,p_wald) using -f argument")
  }
  # annotation: need to have "chr", "locus", "start", "stop","annotation"
  # annot_file = opt$annot_file
  
  if(!is.null(opt$bychrom) & ( isTRUE(opt$bychrom) )){
    BY_CHROMOSOME_PLOTS = TRUE
  }
  
} else {
  # replace with your own values
  trait_name = "Seedling height"
  gwas_result_file = "output/lmm-SeedlHeight.assoc.txt"
  
  # annot_file = "Data/annotations/all.locus_brief_info.7.0"
}





####################
library(readr)
library(qqman)
library(dplyr)

# Load data
gw = read_table2(gwas_result_file)

        if(!("P" %in% names(gw) ) & (! "p_wald" %in% names(gw)) ){
          gw = read_table2(gwfile)
        }
        if("P" %in% names(gw)){
          names(gw)[ names(gw)=="P"]="p_wald"
        }

        if("EMP1" %in% names(gw)){
          names(gw)[ names(gw)=="EMP1"]="p_wald"
        }

        gw$logp = -log10(gw$p_wald)
        if("CHR" %in% names(gw)){
          names(gw)[ names(gw)=="CHR"]="chr"
        }
        if("BP" %in% names(gw)){
          names(gw)[ names(gw)=="BP"]="ps"
        }
        if("SNP" %in% names(gw)){
          names(gw)[ names(gw)=="SNP"]="rs"
        }
gwas=gw

gwas = gwas[ !is.na(gwas$p_wald) ,]

# annotation
#genes = read_tsv(annot_file)
#genes = unique(genes[ genes$is_representative=="Y",
#                      c("chr", "locus", "start", "stop","annotation")])


## QQ plot
#options(bitmaptype="cairo")
png(filename = paste0(trait_name, "-QQplot.png"), type = "cairo-png", 
   width = 4, height = 4, bg = "white", units = "in", res = 250)
   qq(pvector = gwas$p_wald, main=paste0("QQ plot for ", trait_name))
dev.off()


# Manhattan plot
maxY =  max(-log10(gwas$p_wald), na.rm=T)

genomewide= opt$genomewide
if(is.null(genomewide)){
  genomewide = -log10( 0.05/dim(gwas)[[1]])
} else {
  maxY = max(maxY, genomewide)
}
suggestive = opt$suggestive
if(is.null(suggestive)){
  suggestive = 9999
}

#options(bitmaptype='cairo')
png(filename = paste0(trait_name, "-manh-genomewide.png"), 
   type="cairo-png", width = 10, height = 4, bg = "white", 
   units = "in", res = 250)
   qqman::manhattan(gwas, chr = "chr", bp = "ps", snp = "rs", p = "p_wald"   
                 #,annotatePval = 7e-7
                 , main = trait_name
                 , suggestiveline = suggestive
                 , genomewideline = genomewide
                 , ylim =  c(0, 0.6 + maxY)
                 )
dev.off()

#############
# by chromosome

if(BY_CHROMOSOME_PLOTS){
  for(focus_chr in  as.character(unique(gwas$chr)) ){
    
    gwas1 = gwas[ as.character(gwas$chr) %in% focus_chr, ]
    
    xlim = c(0, max(gwas1$ps/1e6) )
    
    png(filename = paste0(trait_name, "-manh-Chr",  focus_chr,".png"), 
        width = 6, height = 4, bg = "white",  units = "in", res = 300)
    qqman::manhattan(gwas1, chr = "chr", 
                     bp = "ps", snp = "rs", 
                     xlim = xlim,
                     p = "p_wald"   
                     #,annotatePval = 7e-7
                     , main = trait_name
                     , suggestiveline = suggestive
                     , genomewideline = genomewide
                     , ylim =  c(0, 0.6 + maxY)
    )
    dev.off()
    
  }
}

#  Possible output:
# gwas[ -log10(gwas$p_wald) > genomewide, c("rs", "chr", "ps")]

# END #

