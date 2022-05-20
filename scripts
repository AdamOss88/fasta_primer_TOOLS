```{r}
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

```
