# FixGeneNames.R
# If gene name starts with number, it would be X"number" because R sucks. Let's fix that.

FixGeneName <- function(gene){
    first.char <- substring(gene, 1, 1)
    if (first.char == "X"){
      second.char <- substring(gene, 2, 2)
      second.char.int <- as.integer(second.char)
      if (!is.na(second.char.int)){
        return(substring(gene, 2, nchar(gene)))
      }
    } else {
      return(gene)
    }
  }