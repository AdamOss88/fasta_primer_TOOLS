```{R}
###function wahtsmissing

whatsmissing <- function(expected,found) {
  
  missing <- NULL
  for (i in 1:length(expected)) {
    if (sum(grepl(expected[i], found)) != 1) {
      missing <- c(missing,expected[i],i)
      i <- i +1
    } else {i <- i +1}
  }
return(missing)  
}

library(seqinr)

###overview GC content and see the length, nr of stop codons
fasta_overview <- function(input) {
  
  n <- NULL    
  m <- NULL    
  
  ##construction
  
  #loop looking for codons
  
  frames_all <- data.frame( stop_F1 = numeric(), stop_F2 = numeric(), stop_F3 = numeric() )
  
  for (cod in 1:length(input)) {
    sequence <- getSequence(input[[cod]])
    
    frames <- c(0,0,0)
    for (frame in 1:3) {
      stop <- 0
      for (m in seq(frame, floor(length(sequence)), 3) ) {
        if (paste(sequence[m:(m+2)], collapse = "" ) == "tag") {
          stop <- stop +1
        } else if (paste(sequence[m:(m+2)], collapse = "" ) == "taa") {
          stop <- stop +1  
        } else if (paste(sequence[m:(m+2)], collapse = "" ) == "tga") {
          stop <- stop +1
        } else {m <- m + 3}
      }
      frames[frame] <- stop
      frame <- frame +1
    }
    frames_all[cod,] <- frames
    
  }  
  
  
  output <- data.frame(annotation = character(), GC = numeric(), length = numeric(),
                       stop_F1 = numeric(), stop_F2 = numeric(), stop_F3 = numeric())
  
  for (n in 1:length(input)) {
    output[n,] <-c(
      getAnnot(input[[n]]),
      GC(getSequence(input[[n]])),
      length(getSequence(input[[n]])),
      frames_all[n,]
    )
  }
  return(output)
}##end of function

hist(as.numeric(geo_genes_overview$length), breaks = 200)
hist(as.numeric(geo_genes_overview$GC), breaks = 200)





###fa_fi - the fasta filter
###
###reads the annotation of fasta file and filter based on a text string
###
###input: fa_fi(input,filter_word)
###input - seqin fasta file object (read.fasta(file = "data/XXXXYYYY.fasta"))
###filter_word - text string or a text strings vector for filtering in .fasta 
###             file name, case sensitive
###output: .fasta files with subset of sequences containing a text string in 
###        the name and a file with remaining sequences not matching any string
###       from filter_word
###usage: fasta_filter(geo_genes,c("Strepto","nost","actino"))
###
fa_fi <- function(input,filter_word) {
  
  require(seqinr)
  positives <- NULL
  rest <- NULL
  n <- NULL
  m <- NULL
  
  if (dir.exists("output_filtering") == F) {
    dir.create("output_filtering")  }
  
  
  for (m in 1:length(filter_word)) {
    
    for (n in 1:length(input)) {
      annotation <- getAnnot(input[n])
      if ( grepl(filter_word[m], annotation, ignore.case = F) ) {
        positives <- append(positives, n)
        n <- n +1
      }
      else {
        n <- n +1
      }
    }
    write.fasta(sequences = input[positives], names = sub(">","",getAnnot(input[positives])), 
                file.out = paste("output_filtering/",filter_word[m],".fasta", sep="")) 
    positives <- NULL
  }
  write.fasta(sequences = input[-rest], names = sub(">","",getAnnot(input[-rest])), 
              file.out = paste("output_filtering/",filter_word[m],"_rest",".fasta", sep=""))
  
  print("possibly done and saved to output_filtering/")
}
###end of function


###filter enzyme names and sizes
###
###
###

##filter names


### filter sizes
### length_filter(input,fraction)
### input - seqinr object , fraction - number 0-1 defining a cutof as a fraction of the mean
### fraction - number 0-1 defining a cutoff as a fraction of the mean (suggested 0,8)

length_filter <- function(input,fraction) {
  require(seqinr)
  #size_treshold <- round(mean(as.numeric(summary(input)[,1]))) * fraction
  size_treshold <- fraction
  i <- NULL
  long_enough <- NULL
  
  for (i in 1:length(input)) {
    if (as.numeric(summary(input)[i,1]) >= size_treshold) {
      long_enough <- append(long_enough , i)
      i <- i+1
    } else {i <- i+1}
    
  }
  
  if (dir.exists("output_filtering") == F) {
    dir.create("output_filtering")  }
  
  write.fasta(sequences = input[long_enough], names = names(input[long_enough]), 
              file.out = paste("output_filtering/",fraction,"_size filter",".fasta", sep="")) 
  output <- read.fasta(file = paste("output_filtering/",fraction,"_size filter",".fasta", sep=""))
  
  print("possibly done and saved to output_filtering/")
  return(output)
  
}###end function



###reverse compliment according to (+/-)strand
### all to (+)
###

All_fasta_53 <- function(input) {
  
  require(seqinr)
  negatives <- NULL
  n <- NULL
  m <- NULL
  
  
  ##get wrongly oriented sequences    
  for (n in 1:length(input)) {
    annotation <- getAnnot(input[[n]])
    if ( grepl("(-)", annotation, fixed = T) ) {
      negatives <- append(negatives, n) 
      n <- n +1
    } else if ( grepl("(+)", annotation, fixed = T) ) {
      n <- n +1
    }
    else {
      print(paste("no information in record ",n))
      n <- n +1
    }
  }
  ##RC them  
  input_seq_rc <- getSequence(input)
  
  for (m in 1:length(negatives)) {
    
    input_seq_rc[[negatives[m]]] <- rev(comp(s2c(c2s(input_seq_rc[[negatives[m]]])))) 
    
  }  
  
  ###save 
  if (dir.exists("output_filtering") == F) {
    dir.create("output_filtering")  }
  
  write.fasta(sequences = input_seq_rc, names = getAnnot(input), 
              file.out = paste("output_filtering/","output_RC",".fasta", sep=""))
  
  output <- read.fasta(file = paste("output_filtering/","output_RC",".fasta", sep=""))  
  return(output)
  
  print("possibly done and saved to output_filtering/")
}###end function

```
