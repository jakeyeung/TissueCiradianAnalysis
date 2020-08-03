# FilterBedFile.R
# February 26 2014

FilterBedFile <- function(genes.in, bedfile.out, bedfile.in){
  source("scripts/functions/ReadListToVector.R")
  source("scripts/functions/FixGeneName.R")
  
  if (missing(bedfile.in)){
    bedfile.in <- "bedfiles/Mm_EPDnew.liftOver.mm10.bed"
  }
  bedfile <- read.table(bedfile.in)  # 4th column contains gene names V4
  genes <- ReadListToVector(genes.in)
  
  # some genes added an 'X' because it R didn't allow rownames with numbers.
  # fix those pesky X0983 ... genenames.
  genes.fixed <- unlist(sapply(genes, FixGeneName))
  
  genes.grep <- paste0(genes.fixed, collapse = "|")
  if (nchar(genes.grep) > 2048){
    warning("Too many genes (>2048 char limit), skipping...")
    return(NA)
  }
  bedfile.filtered <- bedfile[grepl(genes.grep, bedfile$V4), ]
  
  write.table(bedfile.filtered, 
              file = bedfile.out, 
              quote = FALSE, 
              sep = "\t", 
              col.names = FALSE, 
              row.names = FALSE)  
}
