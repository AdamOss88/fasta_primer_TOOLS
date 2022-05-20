---
title: "Primer-utilities"
author: "Adam_Oss"
date: "2/14/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown


```{r degenerate_table}
deg_table <- data.frame(base = c("W","S","M","K","R","Y","B","D","I","H","V","N"),
                        times = c(2,2,2,2,2,2,3,3,3,3,3,4),
                        X1 = c("A","C","A","G","A","C","C","A","A","A","A","A"),
                        X2 = c("T","G","C","T","G","T","G","G","G","C","C","C"),
                        X3 = c(NA,NA,NA,NA,NA,NA,"T","T","T","T","G","T"),
                        X4 = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"G") )
```


```{r function choose degenerate}
choose_degenerate <- function(input){
require(seqinr)
#make a vector of degenerate primers
deg_primers_pos <- NULL
for (i in 1:length(primers)) {
seq <-toupper(getSequence(input[[i]]))
  if (any(deg_table[,1] %in% seq)) {
    deg_primers_pos <- append(deg_primers_pos,i)
    i <- i+1
  } else {
    i <- i+1
  }
rm(seq,i)
}
##get only degenerated
export <- NULL
for (i in 1:length(deg_primers_pos)) {
  export[[i]] <- as.SeqFastadna(input[[deg_primers_pos[i]]], name = sub(">","",getAnnot(input[[deg_primers_pos[i]]])), Annot = getAnnot(input)[[deg_primers_pos[i]]])
}

return(export)

}
###end of function
```

```{r propagate}
degenerateme <- function(input){ 
require(seqinr)
for (pos in 1:length(input)) {
deg_primer <- input[pos]
deg_primer[[1]] <- toupper(deg_primer[[1]])
# has to stay a list to keep one initial primer in one place
for (i in 1:length(deg_primer[[1]])) {
    test_base <- deg_primer[[1]][i]
    if (test_base %in% deg_table[,1]) {
      which <- match(test_base,deg_table[,1]) 
      deg_primer <- append(deg_primer,rep(deg_primer,deg_table[which,2] - 1)) #works until here the nr of copies is ok

#generate a vector with replacements
#papulate
      repl_vect <- sort(rep(as.character(deg_table[which, c(3: (2+deg_table[which,2]))]) ,length(deg_primer) / deg_table[which,2] ))      
      for (j in 1:length(deg_primer)) {deg_primer[[j]][i] <- repl_vect[j]}
      i <- i+1
      rm(test_base, which, repl_vect)
    } else {i <- i+1} 
}
#fix the names as an attribute
names_seq <- NULL
for (i in 1:length(deg_primer)) {
  names_seq[i] <- paste(attr(input[[pos]], which = "name"),"_var_",i, sep = "")  
}

export <- NULL
for (i in 1:length(deg_primer)) {
export[[i]] <- as.SeqFastadna(deg_primer[[i]], name = names_seq[i], Annot = names_seq[i])
}
write.fasta(sequences = getSequence(export), names = getAnnot(export), 
            file.out = paste(attr(input[[pos]], which = "name"),"_variants.fasta",sep = ""))

#rm(export,i,j,test_base,deg_primer)
}
}###end of function
```



