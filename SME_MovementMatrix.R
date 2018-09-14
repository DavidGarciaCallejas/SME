
#' Generate a metacommunity matrix with dispersal and/or foraging behaviours
#'
#' @param N number of communities
#' @param S number of species
#' @param topology "global","linear","random", "lattice". global means distance=1 between any two communities,
#' linear allows movement only connecting adjacent communities, and random creates random connectivity links.
#' "lattice" creates a 2d lattice. 
#' @param ncols if topology is "lattice", number of columns of the lattice. Note that N cannot be and odd number in this case.
#' @param dispersal.coefs "constant","scaling". either constant dispersal rates among species or scaling negatively with niche axis.
#' @param foraging.coefs "constant","scaling". either constant within the foraging range of each sp (constant), or
#' negatively correlated with foraging distance (scaling).
#' @param max.dispersal.distance scalar with the maximum dispersal distance by the species with highest niche value.
#' @param max.foraging.distance scalar with the maximum foraging distance by the species with highest niche value.
#' @param max.dispersal.rate dispersal coefficient for either all species (constant) or that with the lowest niche value (scaling). It is a net dispersal value from
#' a single patch. The effective dispersal to each of N connected patches will be max.dispersal.rate/N
#' @param max.foraging.rate the value of the foraging coefficients relative to the interspecific interaction coefficients. A value of e.g. 0.5 represents
#' that foraging coefficients combined sum up to 0.5 of the interaction coefficient.
#' @param meta.coef.matrix NSxNS matrix with the within-patch interaction coefficients of the metacommunity.
#' @param niche a one dimensional numeric vector representing the niche axis along which maximum foraging distances are scaled.
#'
#' @return list containing two matrices: a NxS matrix with dispersal and/or foraging coefficients, and a NxN matrix with the distances between all pairs of patches
#' @export
#'
#' @examples

# N <- 6
# S <- 3
# topology <- "linear"
# ncols <- 2
# dispersal.coefs <- "constant"
# foraging.coefs <- "constant"
# 
# max.dispersal.distance <- 2
# max.foraging.distance <- 2
# 
# max.dispersal.rate <- 0
# max.foraging.rate <- 0.5
# coef.matrix <- matrix(c(-1.0000000,0.0000000,0.4227851,0.0000000,-1.0000000,0.0000000,-0.1900445,0.0000000,-1.0000000),nrow = 3)
# meta.coef.matrix <- Metacommunity_SME(A = coef.matrix,S = S,N = N,meta.sd.a = 0.1,meta.sd.d = 0.1,type = "top")
# niche <- c(0.2643521,0.6062683,0.9376420)

MovementMatrix_SME <- function(N,S,topology,ncols = 0,dispersal.coefs,foraging.coefs,max.dispersal.distance,max.foraging.distance,max.dispersal.rate,max.foraging.rate,meta.coef.matrix,niche){
  
  # 1 - select the community matrix positions that will be nonzero in each foraged patch
  foraging.template <- matrix(0,S,S)
  
  for(i.row in 1:S){
    for(i.col in 1:S){
      if(meta.coef.matrix[i.row,i.col] < 0){
        if(meta.coef.matrix[i.col,i.row] > 0 & i.col != i.row){
          foraging.template[i.row,i.col] <- -1
          foraging.template[i.col,i.row] <- 1
        }# if antagonism
      }# if negative sign
    }# for i.col
  }# for i.row
  
  # 2 - which patches are connected?
  # fill a separate patch connectivity matrix
  patch.connectivity <- matrix(0,nrow = N,ncol = N)
  
  if(topology == "global"){
    # D = matrix(0,nr=N*S,nc = N*S)
    
    for(i.patch in 0:(N-1)){
      for(j.patch in 0:(N-1)){
        if(i.patch != j.patch){
          # D[((S*i.patch)+1):((S*i.patch)+1+S-1),((S*j.patch)+1):((S*j.patch)+1+S-1)] <- foraging.template
          patch.connectivity[i.patch+1,j.patch+1] <- 1
          patch.connectivity[j.patch+1,i.patch+1] <- 1
        }# if not diagonal block
      }# for j.patch
    }# for i.patch
    
  }else if(topology == "linear"){
    # D = matrix(0,nr=N*S,nc = N*S)
    
    for(i.patch in 0:(N-1)){
      for(j.patch in 0:(N-1)){
        # if adjacent patches or last and first (for creating a torus)
        if(i.patch != j.patch & (abs(i.patch - j.patch) == 1 | abs(i.patch - j.patch) == (N-1))){
          # D[((S*i.patch)+1):((S*i.patch)+1+S-1),((S*j.patch)+1):((S*j.patch)+1+S-1)] <- foraging.template
          patch.connectivity[i.patch+1,j.patch+1] <- 1
          patch.connectivity[j.patch+1,i.patch+1] <- 1
        }# if not diagonal block
      }# for j.patch
    }# for i.patch
    
  }else if(topology == "random"){
    connectivity.prob <- 0.5
    
    # D = matrix(0,nr=N*S,nc = N*S)
    
    for(i.patch in 0:(N-1)){
      for(j.patch in 0:(N-1)){
        if(i.patch != j.patch){
          my.prob <- runif(1,0,1)
          if(my.prob < connectivity.prob){
            # D[((S*i.patch)+1):((S*i.patch)+1+S-1),((S*j.patch)+1):((S*j.patch)+1+S-1)] <- foraging.template
            patch.connectivity[i.patch+1,j.patch+1] <- 1
            patch.connectivity[j.patch+1,i.patch+1] <- 1
          }# if patches connected
        }# if not diagonal block
      }# for j.patch
    }# for i.patch
    
  }else if(topology == "lattice"){
    # D = matrix(0,nr=N*S,nc = N*S)
    nrows <- N/ncols
    
    if(nrows%%1 == 0){
      
      patch.positions <- matrix(1:N,ncol = ncols,nrow = nrows,byrow = T)
      
      # check every position of the patch arrangement
      
      for(i.patch in 1:nrow(patch.positions)){
        for(j.patch in 1:ncol(patch.positions)){
          
          source.patch <- patch.positions[i.patch,j.patch]
          
          if(i.patch == 1){
           
            if(j.patch == 1){
              
              # first row, first column
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              
            }else if(j.patch == ncol(patch.positions)){
              
              # first row, last column
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              
            }else{
              
              # first row, middle column
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              
            }
             
          }else if(i.patch == nrow(patch.positions)){
            
            if(j.patch == 1){
              
              # last row, first column
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              
            }else if(j.patch == ncol(patch.positions)){
              
              # last row, last column
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              
            }else{
              
              # last row, middle column
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              
            }
            
          }else{
            
            if(j.patch == 1){
              
              # middle row, first column
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              
            }else if(j.patch == ncol(patch.positions)){
              
              # middle row, last column
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              
            }else{
              
              # middle row, middle column
              patch.connectivity[source.patch,patch.positions[i.patch+1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch-1,j.patch]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch+1]] <- 1
              patch.connectivity[source.patch,patch.positions[i.patch,j.patch-1]] <- 1
              
              patch.connectivity[patch.positions[i.patch+1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch-1,j.patch],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch+1],source.patch] <- 1
              patch.connectivity[patch.positions[i.patch,j.patch-1],source.patch] <- 1
              
            }
            
          }# if-else i.patch ==1
        } # for j.patch
      } #for i.patch
      
    }else{
      stop("provide appropriate N values for the lattice topology")
    }# if correct N
  }# topology
  
  # 3 - generate matrix of distances between patches
  # igraph structure from adjacency matrix
  patch.graph <- igraph::graph_from_adjacency_matrix(adjmatrix = patch.connectivity,mode = "undirected",diag = FALSE)
  # calculate distances between patches
  patch.distances <- igraph::distances(graph = patch.graph,mode = "out")
  
  # 4 - foraging and dispersal niches
  
  # remap the niches to the range 0-d
  foraging.niche <- c(0,niche,1)
  foraging.niche <- foraging.niche - min(foraging.niche)
  foraging.niche <- foraging.niche/max(foraging.niche)
  foraging.niche <- foraging.niche * max.foraging.rate
  foraging.niche <- foraging.niche[2:(length(foraging.niche)-1)]
  
  dispersal.niche <- c(0,niche,1)
  dispersal.niche <- dispersal.niche - min(dispersal.niche)
  dispersal.niche <- dispersal.niche/max(dispersal.niche)
  dispersal.niche <- dispersal.niche * max.dispersal.rate
  dispersal.niche <- dispersal.niche[2:(length(dispersal.niche)-1)]
  
  # create metacommunity matrix
  D = matrix(0,nr=N*S,nc = N*S)
  
  # 1. maximum distances reachable by any species
  
  if(max.dispersal.rate > 0){
    
  sp.max.distance <- (max.dispersal.distance+0.1)/max(dispersal.niche)*dispersal.niche
  # matrix containing the dispersal coefficient of each sp for each distance
  dispersal.coef.data <- matrix(0,nrow=S,ncol=max(patch.distances))
  
  if(dispersal.coefs == "constant"){
    # every species disperses at least the nearest patch
    dispersal.coef.data[,1] <- max.dispersal.rate
    
    # the rest of the coefficients are constant within the reachable distances
    for(i.row in 1:nrow(dispersal.coef.data)){
      if(sp.max.distance[i.row]>2){
        dispersal.coef.data[i.row,2:floor(sp.max.distance[i.row])] <- max.dispersal.rate
      }# if reaches other patches
    }# for each row
    
    # divide dispersal across reachable patches
    if(ncol(dispersal.coef.data)>1){
      dispersal.coef.data <- t(apply(X = dispersal.coef.data,MARGIN = 1,FUN = function(x) x/sum(x>0)))
    }
    
  }else if(dispersal.coefs == "scaling"){
    
    # this formula is a line equation expanded to 1) find m and 2) find b, and these plugged into y=mx*b to get the values of dispersal(y) for a certain distance(x)
    for(i.sp in 1:nrow(dispersal.coef.data)){
      my.sp.distances <- ((max.dispersal.rate/(1-sp.max.distance[i.sp]))*1:floor(sp.max.distance[i.sp])) + (max.dispersal.rate*(-sp.max.distance[i.sp]))/(1-sp.max.distance[i.sp])
      if(length(my.sp.distances)<ncol(dispersal.coef.data)){
        my.sp.distances[(length(my.sp.distances)+1):ncol(dispersal.coef.data)] <- 0
      }
      dispersal.coef.data[i.sp,1:ncol(dispersal.coef.data)] <- my.sp.distances[1:ncol(dispersal.coef.data)] 
    }# for i.sp
    dispersal.coef.data[dispersal.coef.data<0] <- 0
    
    # divide dispersal effort across reachable patches
    for(i.row in 1:nrow(dispersal.coef.data)){
      if(sum(dispersal.coef.data[i.row,]>0)>1){
        dispersal.coef.data[i.row,] <- dispersal.coef.data[i.row,] * (max.dispersal.rate/sum(dispersal.coef.data[i.row,]))
      }# if more than one patch
    }# for i.row
    
  }# if-else dispersal is constant or scaling
  
  # 5 - assign dispersal coefficients
  # maintaining mass balance, i.e. the sum of the dispersal to all connected patches i.e. sum(d/N-1)
  # is equivalent to the loss of the source patch, i.e. -d
  # -d is the value of the diagonal patch, then.

    # update matrix D
    for(i.row in 1:nrow(D)){
      for(i.col in 1:ncol(D)){
        
        # if diagonal element and not main diagonal,
        # it is a dispersal coefficient
        if(i.row %% S == i.col %% S & i.row != i.col){
          
          # which sp?
          my.sp <- i.row %% S
          if(my.sp == 0){my.sp <- S}
          
          # source and dest patch
          source.patch <- ceiling(i.col/S) 
          dest.patch <- ceiling(i.row/S)
          
          # double check
          if(source.patch != dest.patch){
            
            # distance among them
            my.distance <- patch.distances[source.patch,dest.patch]
            
            # the value is the one from dispersal.coef.data (the net effort for a given distance)
            # divided by the number of patches at that distance
            
            # how many patches at this distance?
            patches.at.dist <- table(patch.distances[source.patch,])
            patches.at.dist <- as.integer(patches.at.dist[which(names(patches.at.dist) == my.distance)])
            
            # divide the effort from foraging.coef.data among all patches within a given distance
            dispersal.coef.weight <- dispersal.coef.data[my.sp,my.distance]/patches.at.dist
            
            D[i.row,i.col] <- dispersal.coef.weight
            
          }# if different patch
        }# if dispersal cell
        
      }# for i.col
    }# for i.row
  
  # 6 - diagonal elements reflect the loss of population due to dispersal
  # (see section on dispersal coefficients)
  diag(D) <- -max.dispersal.rate
  
  }# if max.dispersal.rate > 0
  
  # 7 - assign foraging coefficients
  
  # 7.1 - calculate the foraging distance and coefficients for every species
  # considering that species differences in their niche axis (body size)
  # are related to their foraging distance (Jetz et al. 2003)
  
  # assume that species with the highest niche value forages up to the maximum distance (and a bit more,0.4, to ensure a value >0),
  # and a linear relationship for estimating the other species' distances
  # the linear relationship crossess (0,0)
  # I also assume that all species forage to the nearest patch with max.foraging.rate coefficient
  # regardless of their niche value
  
  # note also that with foraging there is no dispersal associated, and no mass balance constraints
  # ASIDE FROM -POTENTIALLY- decreasing pressure over source patch
  # in the first version coded, max.foraging.rate is relative to the interspecific interaction coefficients
  # so that e.g. max.foraging.rate = 1 means that the overall foraging is equal to the a_ij coefficient of the community
  
  if(max.foraging.rate > 0){
    
      sp.max.distance <- (max.foraging.distance+0.1)/max(foraging.niche)*foraging.niche
      # matrix containing the foraging coefficient of each sp for each distance
      foraging.coef.data <- matrix(0,nrow=S,ncol=max(patch.distances))
      
      if(foraging.coefs == "constant"){
        # every species forages at least the nearest patch
        foraging.coef.data[,1] <- max.foraging.rate
        
        # the rest of the coefficients are constant within the reachable distances
        for(i.row in 1:nrow(foraging.coef.data)){
          if(sp.max.distance[i.row]>2){
            foraging.coef.data[i.row,2:floor(sp.max.distance[i.row])] <- max.foraging.rate
          }# if reaches other patches
        }# for each row
        
        # divide foraging effort across reachable patches
        if(ncol(foraging.coef.data)>1){
          foraging.coef.data <- t(apply(X = foraging.coef.data,MARGIN = 1,FUN = function(x) x/sum(x>0)))
        }
        
      }else if(foraging.coefs == "scaling"){
        
        # this formula is a line equation expanded to 1) find m and 2) find b, and these plugged into y=mx*b to get the values of foraging(y) for a certain distance(x)
        for(i.sp in 1:nrow(foraging.coef.data)){
          my.sp.distances <- ((max.foraging.rate/(1-sp.max.distance[i.sp]))*1:floor(sp.max.distance[i.sp])) + (max.foraging.rate*(-sp.max.distance[i.sp]))/(1-sp.max.distance[i.sp])
          if(length(my.sp.distances)<ncol(foraging.coef.data)){
            my.sp.distances[(length(my.sp.distances)+1):ncol(foraging.coef.data)] <- 0
          }
          foraging.coef.data[i.sp,1:ncol(foraging.coef.data)] <- my.sp.distances[1:ncol(foraging.coef.data)] 
        }# for i.sp
        foraging.coef.data[foraging.coef.data<0] <- 0
        
        # divide foraging effort across reachable patches
        for(i.row in 1:nrow(foraging.coef.data)){
          if(sum(foraging.coef.data[i.row,]>0)>1){
            foraging.coef.data[i.row,] <- foraging.coef.data[i.row,] * (max.foraging.rate/sum(foraging.coef.data[i.row,]))
          }# if more than one patch
        }# for i.row
      }# if-else coefs

    
    # update matrix D
    
    # see that these statements work:
    # for(i.col in 1:ncol(D)) print(paste("col:",i.col,"-patch:",ceiling(i.col/S),"-sp:",ifelse(i.col %% S == 0,S,i.col%%S)))
    
    for(i.row in 1:nrow(D)){
      for(i.col in 1:ncol(D)){
        
        forager.sp <- ifelse(i.col %% S == 0,S,i.col %% S)
        prey.sp <- ifelse(i.row %% S == 0,S,i.row %% S)
        
        # is there predation?
        if(foraging.template[prey.sp,forager.sp] == -1){
          
          # if so, check distance between patches
          source.patch <- ceiling(i.col/S) 
          dest.patch <- ceiling(i.row/S)
          
          if(source.patch != dest.patch){
            my.distance <- patch.distances[source.patch,dest.patch]
            
            # how many patches at this distance?
            patches.at.dist <- table(patch.distances[source.patch,])
            patches.at.dist <- as.integer(patches.at.dist[which(names(patches.at.dist) == my.distance)])
            
            # divide the effort from foraging.coef.data among all patches within a given distance
            coef.weight <- foraging.coef.data[forager.sp,my.distance]/patches.at.dist
            
            # final value is the specific weight multiplied by the original coefficient
            D[i.row,i.col] <- coef.weight * meta.coef.matrix[(source.patch-1)*S+prey.sp,(source.patch-1)*S+forager.sp]
            D[i.col,i.row] <- coef.weight * meta.coef.matrix[(source.patch-1)*S+forager.sp,(source.patch-1)*S+prey.sp]
            
          }# if different patch
        }# if predation
        
      }# for i.col
    }# for i.row

    # substract the foraging effort from the focal patch coefficient. In other words,
    # if a species forages in other patches, its effect on local prey will be smaller
    
    for(i.patch in 1:N){
      my.row.positions <- S*(i.patch-1) + 1:S
      
      my.positions <- expand.grid(my.row.positions,my.row.positions)
      my.positions <- my.positions[my.positions[,1] != my.positions[,2],]
      
      # D[my.row.positions,my.col.positions] <- -max.foraging.rate * meta.coef.matrix[my.row.positions,my.col.positions]
      D[cbind(my.positions[,1],my.positions[,2])] <- max.foraging.rate * meta.coef.matrix[cbind(my.positions[,1],my.positions[,2])]
      
    }# for i.patch
  }# if max.foraging.rate > 0
  
  return(list(D,patch.distances))
}

