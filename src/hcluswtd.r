dissim <- function(a, wt) {
# Inputs.   a: matrix, for which we want distances on rows,
#           wt: masses of each row.
# Returns.  matrix of dims. nrow(a) x nrow(a) with wtd. sqd. Eucl. distances.
# FM, 2003/11/16

   n <- nrow(a)
   m <- ncol(a)
   adiss <- matrix(0, n, n)

   for (i1 in 2:n) {
       adiss[i1,i1] <- 0.0
       for (i2 in 1:(i1-1)) {
           adiss[i1,i2] <- 0.0
           for (j in 1:m) {
               # We use the squared Euclidean distance, weighted.
               adiss[i1,i2] <- adiss[i1,i2] + (wt[i1]*wt[i2])/(wt[i1]+wt[i2]) *
                    (a[i1,j]-a[i2,j])^2
           }
           adiss[i2,i1] <- adiss[i1,i2]
       }
   }
   adiss
}


getnns <- function(diss, flag) {
# Inputs.  diss: full distance matrix.
#          flag: "live" rows indicated by 1 are to be processed.
# Returns. List of: nn, nndiss.
#          nn:   list of nearest neighbor of each row.
#          nndiss: nearest neigbbor distance of each row.
# FM, 2003/11/16

   nn <- rep(0, nrow(diss))
   nndiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   if (nrow(diss) != ncol(diss)) stop("Invalid input first parameter.")
   if (nrow(diss) != length(flag)) stop("Invalid inputs 1st/2nd parameters.")
 # if (nrow(diss) != length(nn)) stop("Invalid inputs 1st/3rd parameters.")
 # if (nrow(diss) != length(nndiss)) stop("Invalid inputs 1st/4th parameters.")

   for (i1 in 1:nrow(diss)) {
       if (flag[i1] == 1) {
          minobs <- -1
          mindis <- MAXVAL
          for (i2 in 1:ncol(diss)) {
              if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
                 mindis <- diss[i1,i2]
                 minobs <- i2
              }
          }
          nn[i1] <- minobs
          nndiss[i1] <- mindis
       }
   }
   list(nn = nn, nndiss = nndiss)
}

hierclust <- function(a, wt) {

   MAXVAL <- 1.0e12

   n <- nrow(a)                               
   diss <- dissim(a, wt)                      # call to function dissim
   flag <- rep(1, n)                          # active/dead indicator
   a <- rep(0, n-1)                           # left subnode on clustering
   b <- rep(0, n-1)                           # right subnode on clustering
   ia <- rep(0, n-1)                          # R-compatible version of a
   ib <- rep(0, n-1)                          # R-compatible version of b
   lev <- rep(0, n-1)                         # level or criterion values
   card <- rep(1, n)                          # cardinalities
   mass <- wt
   order <- rep(0, n)                         # R-compatible order for plotting

   nnsnnsdiss <- getnns(diss, flag)           # call to function getnns
   clusmat <- matrix(0, n, n)                 # cluster memberships
   for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition

   for (ncl in (n-1):1) {                      # main loop 
       # check for agglomerable pair
       minobs <- -1;  
       mindis <- MAXVAL;
       for (i in 1:n) {
           if (flag[i] == 1) {
              if (nnsnnsdiss$nndiss[i] < mindis) {
                  mindis <- nnsnnsdiss$nndiss[i]
                  minobs <- i
              }
           }
       }
       # find agglomerands clus1 and clus2, with former < latter
       if (minobs < nnsnnsdiss$nn[minobs]) {
          clus1 <- minobs
          clus2 <- nnsnnsdiss$nn[minobs]
       }
       if (minobs > nnsnnsdiss$nn[minobs]) {
          clus2 <- minobs
          clus1 <- nnsnnsdiss$nn[minobs]
       }
       # So, agglomeration of pair clus1 < clus2 defines cluster ncl

       #------------------------------------ Block for subnode labels 
       a[ncl] <- clus1                       # aine, or left child node
       b[ncl] <- clus2                       # benjamin, or right child node
       # Now build up ia, ib as version of a, b which is R-compliant
       if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
       if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
       if (card[clus1] > 1) {                # left child is non-singleton
          lastind <- 0
          for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
              if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
          }
          ia[ncl] <- n - lastind             # label of non-singleton
       }
       if (card[clus2] > 1) {                # right child is non-singleton
          lastind <- 0
          for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
              if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
          }
          ib[ncl] <- n - lastind             # label of non-singleton
       }
       if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
          left <- min(ia[ncl],ib[ncl])
          right <- max(ia[ncl],ib[ncl])
          ia[ncl] <- left                    # Just get left < right
          ib[ncl] <- right
       }
       #--------------------------------------------------------------------

       lev[ncl] <- mindis
       for (i in 1:n) {
           clusmat[i,ncl] <- clusmat[i,ncl+1]
           if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
       }
       # Next we need to update diss array
       for (i in 1:n) {
           if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
              diss[clus1,i] <- 
      ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,i] +
      ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus2,i] -
      (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,clus2] 
              diss[i,clus1] <- diss[clus1,i]
           }
       }
       mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
       card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
       # Cluster label clus2 is knocked out; following not nec. but no harm
       flag[clus2] <- 0
       nnsnnsdiss$nndiss[clus2] <- MAXVAL
       mass[clus2] <- 0.0
       for (i in 1:n) {
           diss[clus2,i] <- MAXVAL
           diss[i,clus2] <- diss[clus2,i]
       }
       # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
       # i.e. nearest neighbors and the nearest neigh. dissimilarity
       nnsnnsdiss <- getnns(diss, flag)
   }

   temp <- cbind(a,b)
   merge2 <- temp[nrow(temp):1, ]
   temp <- cbind(ia,ib)
   merge <- temp[nrow(temp):1,]
   dimnames(merge) <- NULL
   # merge is R-compliant; later suppress merge2

   #-------------------------------- Build R-compatible order from ia, ib
   orderlist <- c(merge[n-1,1], merge[n-1,2])
   norderlist <- 2
   for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
       for (i2 in 1:norderlist) {       # Scan orderlist
           if (orderlist[i2] > 0) {     # Non-singleton to be expanded
              tobeexp <- orderlist[i2]
              if (i2 == 1) {
                 orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                                orderlist[2:norderlist])
              }
              if (i2 == norderlist) {
                 orderlist <- c(orderlist[1:(norderlist-1)],
                                merge[tobeexp,1],merge[tobeexp,2])
              }
              if (i2 > 1 && i2 < norderlist) {
                 orderlist <- c(orderlist[1:(i2-1)], 
                                merge[tobeexp,1],merge[tobeexp,2],
                                orderlist[(i2+1):norderlist])
              }
              norderlist <- length(orderlist)
           }
        }
   }
   orderlist <- (-orderlist)
   class(orderlist) <- "integer"
   
   xcall <- "hierclust(a,wt)"
   class(xcall) <- "call"
   #clusmat=clusmat
   #labels=as.character(1:n)

   retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
         labels=dimnames(a)[[1]],method="minvar",call=xcall,
         dist.method="euclidean-factor")
   retlist <- list(merge=merge,height=lev[(n-1):1],order=orderlist)
   class(retlist) <- "hclust"
   retlist
}
