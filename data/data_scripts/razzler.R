

# take a ped
# return a tped

razzler <- function(a,b) {
  l <- vector("list", length(a)+length(b))

  j <- 1
  for (i in 1:length(a)) {
    l[[j]] <- a[i]
    l[[j+1]] <- b[i]
    j <- j+2
  }
  return(unlist(l))
}


makeATped <- function(ped, fileprefix) {

    n <- ncol(ped)
    d <- ped[,7:n] # just the genotypes
    m <- ncol(d)
    
    tped <- matrix(nrow=(m/2), ncol=(2*nrow(ped)), data=0)

    for (i in 1:(m/2)) {
        tped[i,] <- razzler(d[,(2*i-1)], d[,2*i])
    }

    fam <- ped[,1:6]
    fone <- paste(fileprefix, "tped", sep="")
    ftwo <- paste(fileprefix, "fam", sep="")
    
    write.table(tped, file=fone, quote=F, row.names=F, col.names=F, sep="\t")
    write.table(fam,  file=ftwo, quote=F, row.names=F, col.names=F, sep="\t")

}
