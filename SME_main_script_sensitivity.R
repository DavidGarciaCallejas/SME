##########################
# main script, but saving the results in order to perform sensitivity analyses
##########################
# it is not very user-friendly, it should be run separately for each dispersal and/or foraging value
# this was done to have an arbitrary number of r sessions in parallel

# write the results to a file?
store.results <- TRUE

# write only the results for these values
dispersal.value <- 0.75
foraging.value <- 0.75

# replicates
replicates <- 50

# verbose output?
verbose <- TRUE

# how many steps from min to max dispersal/foraging?
nstep <- 100

# load parameters 
#######################
#######################
# set of simulations
# 1 - global, constant dispersal
# 2 - lattice, constant dispersal
# 3 - random, constant dispersal
# 4 - lattice, constant foraging
# 5 - lattice, scaling foraging
# 6 - as 2, with only antagonism, c = 0.2
# 7 - as 4, with only antagonism, c = 0.2
# 8 - as 5, with only antagonism, c = 0.2
#######################
# 9 - lattice, constant dispersal,constant foraging
# 10 - lattice, scaling dispersal, scaling foraging
# 11 - as 9, with only antagonism, c = 0.2
# 12 - as 10, with only antagonism, c = 0.2
#######################

ID <- 11

pars = read.table("./data/parameters_set.csv",header = T,sep = ";",dec = ".",stringsAsFactors = F)

# varying dispersal and/or foraging rates
dispersal.rate.min <- 1e-4
foraging.rate.min <- 1e-6

###########
S <- pars$S[ID]	
connectance.type <- pars$connectance.type[ID]
connectances <- c(pars$amensalism.connectance[ID],
                  pars$antagonism.connectance[ID],
                  pars$commensalism.connectance[ID],
                  pars$competition.connectance[ID],
                  pars$mutualism.connectance[ID])
interaction.structures <- c(pars$amensalism[ID],
                            pars$antagonism[ID],
                            pars$commensalism[ID],
                            pars$competition[ID],
                            pars$mutualism[ID])

mean.a <- pars$mean.a[ID]
sd.a <- pars$sd.a[ID]
mean.d <- pars$mean.d[ID]
sd.d <- pars$sd.d[ID]
N <- pars$N[ID]
meta.sd.a <- pars$meta.sd.a[ID]
meta.sd.d <- pars$meta.sd.d[ID]
meta.type <- pars$type[ID]
topology <- pars$topology[ID]#"global","linear","random","lattice"
topology.cols <- 5
dispersal.coefs <- pars$dispersal.coefs[ID]#"constant","scaling"
foraging.coefs <- pars$foraging.coefs[ID]
max.dispersal.rate <- ifelse(pars$max.dispersal.rate[ID] > 0, dispersal.value,0)
max.foraging.rate <- ifelse(pars$max.foraging.rate[ID] > 0, foraging.value,0)

###########
max.dispersal.distance <- 2
max.foraging.distance <- 2
###########

dispersal.param <- ifelse(max.dispersal.rate > 0,"intermediate_disp","no_disp")
if(max.dispersal.rate > 0){
  
  if(dispersal.value > 0.7){
    dispersal.param <- "high_disp"
  }else if(dispersal.value < 0.3){
    dispersal.param <- "low_disp"
  }
}

foraging.param <-  ifelse(max.foraging.rate > 0,"intermediate_for","no_for")
if(max.foraging.rate > 0){
  
  if(foraging.value > 0.7){
    foraging.param <- "high_for"
  }else if(foraging.value < 0.3){
    foraging.param <- "low_for"
  }
}

ID <- paste(ID,"_",dispersal.param,"_",foraging.param,sep = "")

##########################
##########################
# auxiliary functions
get.distances <- function(x,y){
  patch.distances[x,y]
}

get.path.lengths <- function(x,y){
  path.lengths[x,y]
}

##########################
##########################
# results dataframe
results.data <- data.frame(simulation.ID = rep(ID,1e6),
                           replicate = 0,
                           dispersal.rate = 0,
                           foraging.rate = 0,
                           affected.sp = 0,
                           trigger.sp = 0,
                           affected.sp.patch = 0,
                           trigger.sp.patch = 0,
                           spatial.distance = 0,
                           path.length = 0,
                           lambda.max = 0,
                           direct.effect = 0,
                           net.effect = 0,stringsAsFactors = FALSE)

trophic.info <- NULL

# couple of auxiliary variables
results.pos <- 1
written <- 1

##########################
##########################
for(i.replicate in 1:replicates){
  
  niche <- sort(runif(S,0,1))
  
  ##########################
  ##########################
  
  # Find the starting coefficients matrix
  # and make sure there are no isolated components
  
  test=numeric(S*N)+1
  num.components <- 2
  while(sum(test)>0 | num.components > 1) {
    test=numeric(S*N)
    sign.matrix <- SignMatrix_SME(S = S,
                                  connectances = connectances,
                                  amensalism = interaction.structures[1],
                                  antagonism = interaction.structures[2],
                                  commensalism = interaction.structures[3],
                                  competition = interaction.structures[4],
                                  mutualism = interaction.structures[5],
                                  connectance.type = connectance.type,
                                  niche = niche)
    
    A <- InteractionCoefs_SME(sign.matrix = sign.matrix,mean.a = mean.a,sd.a = sd.a,mean.d = mean.d,sd.d = sd.d,symmetric = TRUE)
    metaA <- Metacommunity_SME(A = A,S = S,N = N,meta.sd.a = meta.sd.a,meta.sd.d = meta.sd.d,type = meta.type)
    SS = eq_local(metaA,S,N,nstep=10)
    test[SS==0]=1
    
    if(sum(test)==0){
      
      # calculate trophic positions
      antagonism.matrix <- matrix(0,S,S)
      
      for(i.row in 1:nrow(A)){
        for(i.col in 1:ncol(A)){
          if(A[i.row,i.col] < 0){
            if(A[i.col,i.row] > 0 & i.col != i.row){
              antagonism.matrix[i.row,i.col] <- A[i.row,i.col]
              antagonism.matrix[i.col,i.row] <- A[i.col,i.row]
            }# if antagonism
          }# if negative sign
        }# for i.col
      }# for i.row
      
      # igraph connected components
      graph.A <- igraph::graph_from_adjacency_matrix(adjmatrix = A,mode = "undirected",weighted = "1",diag = FALSE)
      num.components <- components(graph.A)$no
      
      if(num.components == 1){
        
        # calculate trophic positions
        # for that, the easiest way is to accomodate food web to the cheddar package
        nodes <- data.frame(node = paste("sp.",1:S,sep=""),stringsAsFactors = F)
        colnames(antagonism.matrix) <- nodes$node
        rownames(antagonism.matrix) <- nodes$node
        expanded.links <- as.data.frame.table(antagonism.matrix,stringsAsFactors = FALSE)
        expanded.links <- subset(expanded.links, Freq < 0)
        names(expanded.links) <- c("resource","consumer","coef")
        my.food.web <- Community(nodes = nodes,properties=list(title='foo'),trophic.links = expanded.links)
        trophic.positions <- data.frame(species = 1:S, PreyAvg.TL = PreyAveragedTrophicLevel(my.food.web))
        
      }# if fully connected
    }# if valid abundances
  }# while zero abundances or isolated components
  
  ##########################
  ##########################
  # variable dispersal and/or foraging rates
  # note that increment is not linear, so further sampling from these rates should not be linear
  
  movement.rates <- data.frame(dispersal.rates = rep(0,nstep),foraging.rates = rep(0,nstep))
  if(max.dispersal.rate != 0){
    movement.rates$dispersal.rates <- 10^seq(log(dispersal.rate.min,10),log(max.dispersal.rate,10),(log(max.dispersal.rate,10)-log(dispersal.rate.min,10))/(nstep-1))
    # just to be sure it's exactly that value
    movement.rates$dispersal.rates[nrow(movement.rates)] <- max.dispersal.rate
    
  }
  if(max.foraging.rate != 0){
    movement.rates$foraging.rates <- 10^seq(log(foraging.rate.min,10),log(max.foraging.rate,10),(log(max.foraging.rate,10)-log(foraging.rate.min,10))/(nstep-1))
    # just to be sure it's exactly that value
    movement.rates$foraging.rates[nrow(movement.rates)] <- max.foraging.rate
    
  }
  
  closest.dispersal.value <- movement.rates$dispersal.rates[vapply(dispersal.value,function(x) which.min(abs(movement.rates$dispersal.rates-x)),1)]
  closest.foraging.value <- movement.rates$foraging.rates[vapply(foraging.value,function(x) which.min(abs(movement.rates$foraging.rates-x)),1)]
  
  ##########################
  ##########################
  # solve the system and store results
  for(i.step in 1:nrow(movement.rates)) {
    
    # generate connectivity matrices
    dispersal.matrices <- MovementMatrix_SME(N = N,
                                             S = S,
                                             topology = topology,
                                             ncols = topology.cols,
                                             dispersal.coefs = dispersal.coefs,
                                             foraging.coefs = foraging.coefs,
                                             max.dispersal.distance = max.dispersal.distance,
                                             max.foraging.distance = max.foraging.distance,
                                             max.dispersal.rate = movement.rates$dispersal.rates[i.step],
                                             max.foraging.rate = movement.rates$foraging.rates[i.step],
                                             meta.coef.matrix = metaA,
                                             niche = niche)
    D <- dispersal.matrices[[1]]
    patch.distances <- dispersal.matrices[[2]]
    
    # path lengths
    graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = metaA+D,mode = "undirected",weighted = "1",diag = FALSE)
    path.lengths <- distances(graph = graph.D,algorithm = "unweighted")
    
    # Steady state solution
    SS = stode(y = SS,time=0,func=model_stode,parms=NULL,A=metaA,D=D, positive = TRUE)[[1]]
    
    # Jacobian
    J = jacobian.full(y=SS,func=model_J,A=metaA,D=D)
    
    # Local stability analysis
    lmax = max(as.double(eigen(J)$values))	
    
    # negative of the inverse jacobian: net effects matrix
    invJ <- -ginv(J)
    
    ################
    ################
    # store results
    
    if((movement.rates$dispersal.rate[i.step] == closest.dispersal.value & closest.dispersal.value != 0) | 
       (movement.rates$foraging.rate[i.step] == closest.foraging.value & closest.foraging.value != 0)){
      
      if(store.results){
        
        # trophic info data
        # don't bother memory management for this dataframe. 
        # it's small enough
        trophic.positions$replicate <- i.replicate
        trophic.positions$niche.axis <- niche
        
        trophic.info <- bind_rows(trophic.info,trophic.positions)
        
        # expand the jacobian
        colnames(J) <- 1:(N*S)
        rownames(J) <- 1:(N*S)
        expanded.jacobian <- as.data.frame.table(J,stringsAsFactors = FALSE)
        names(expanded.jacobian) <- c("affected.sp","trigger.sp","direct.effect")
        expanded.jacobian$trigger.sp <- as.integer(expanded.jacobian$trigger.sp)
        expanded.jacobian$affected.sp <- as.integer(expanded.jacobian$affected.sp)
        
        expanded.jacobian$trigger.sp.patch <- ceiling(expanded.jacobian$trigger.sp/S) 
        expanded.jacobian$affected.sp.patch <- ceiling(expanded.jacobian$affected.sp/S) 
        expanded.jacobian$trigger.sp <- ifelse(expanded.jacobian$trigger.sp %% S == 0,S,expanded.jacobian$trigger.sp %% S)
        expanded.jacobian$affected.sp <- ifelse(expanded.jacobian$affected.sp %% S == 0,S,expanded.jacobian$affected.sp %% S)
        
        # expand the inverse jacobian
        colnames(invJ) <- 1:(N*S)
        rownames(invJ) <- 1:(N*S)
        expanded.invJ <- as.data.frame.table(invJ,stringsAsFactors = FALSE)
        names(expanded.invJ) <- c("affected.sp","trigger.sp","net.effect")
        expanded.invJ$trigger.sp <- as.integer(expanded.invJ$trigger.sp)
        expanded.invJ$affected.sp <- as.integer(expanded.invJ$affected.sp)
        
        expanded.invJ$trigger.sp.patch <- ceiling(expanded.invJ$trigger.sp/S) 
        expanded.invJ$affected.sp.patch <- ceiling(expanded.invJ$affected.sp/S) 
        expanded.invJ$trigger.sp <- ifelse(expanded.invJ$trigger.sp %% S == 0,S,expanded.invJ$trigger.sp %% S)
        expanded.invJ$affected.sp <- ifelse(expanded.invJ$affected.sp %% S == 0,S,expanded.invJ$affected.sp %% S)
        
        expanded.invJ <- dplyr::left_join(expanded.invJ,expanded.jacobian, by=c("affected.sp", "trigger.sp", "trigger.sp.patch", "affected.sp.patch"))
        
        expanded.invJ$spatial.distance <- mapply(FUN = get.distances,expanded.invJ$trigger.sp.patch,expanded.invJ$affected.sp.patch)
        expanded.invJ$path.length <- mapply(FUN = get.path.lengths,
                                            (expanded.invJ$trigger.sp.patch-1)*S+expanded.invJ$trigger.sp,
                                            (expanded.invJ$affected.sp.patch-1)*S+expanded.invJ$affected.sp)
        
        expanded.invJ$dispersal.rate <- movement.rates$dispersal.rates[i.step]
        expanded.invJ$foraging.rate <- movement.rates$foraging.rates[i.step]
        
        expanded.invJ$lambda.max <- lmax
        expanded.invJ$simulation.ID <- ID
        expanded.invJ$replicate <- i.replicate
        
        # trim results
        expanded.invJ <- expanded.invJ[,c("simulation.ID",#"N","S","topology",
                                          "replicate",
                                          "dispersal.rate",
                                          "foraging.rate",
                                          "trigger.sp",
                                          "affected.sp",
                                          "trigger.sp.patch",
                                          "affected.sp.patch",
                                          "spatial.distance",
                                          "path.length",
                                          "lambda.max",
                                          "direct.effect","net.effect")]
        
        # add them to the preallocated big dataframe
        results.data[results.pos:(results.pos+nrow(expanded.invJ)-1),] <- expanded.invJ
        results.pos <- results.pos + nrow(expanded.invJ)
        
        # this is a cheap hack not to overflow the results.data dataframe
        if(results.pos > 990000 | i.replicate == replicates){
          readr::write_delim(x = results.data[1:results.pos,],path = paste("./results/results_replicates_ID",ID,"_file",sprintf("%02d",written),".csv",sep=""),delim = ";",append = F)
          written <- written + 1
          results.pos <- 1
        }
      }# if store.results
      
      if(verbose){
        trimmed.invJ <- invJ
        
        # set 0 to diagonal blocks
        diagonals::fatdiag(trimmed.invJ,S) <- 0
        # set 0 to dispersal rates
        for(i in 1:nrow(trimmed.invJ)){
          for(j in 1:ncol(trimmed.invJ)){
            if(i%%S == j%%S){
              trimmed.invJ[i,j] <- 0
            }
          }# for j
        }# for i
        
        cat("replicate = ", i.replicate,
            "- dispersal:",round(movement.rates$dispersal.rates[i.step],5),
            "- foraging:",round(movement.rates$foraging.rates[i.step],5),
            "- lmax = ",round(lmax,5),
            "- mean -:",round(mean(trimmed.invJ[trimmed.invJ<0]),5),
            "- mean +:",round(mean(trimmed.invJ[trimmed.invJ>0]),5),'\n')
      }
    }# if dispersal.value or foraging.value  
    
  }# for i.step	
  
}# for i.replicate

if(store.results){
  print("storing results...")
  readr::write_delim(x = trophic.info,path = paste("./results/species_data_replicates_ID",ID,".csv",sep=""),delim = ";",append = F)
  
  ## rearrange results in one single file
  results.list <- list()
  for(i.ID in 1:length(ID)){
    num.files <- list.files(path = "./results",pattern = paste("results_replicates_ID",ID[i.ID],"_file",sep=""))
    for(i.file in 1:length(num.files)){
      results.list[[i.file]] <- readr::read_delim(paste("./results/",num.files[i.file],sep=""),delim = ";",col_types = cols())
    }
  }
  full.data <- bind_rows(results.list)
  full.data <- droplevels(subset(full.data,replicate != 0))
  
  readr::write_delim(x = full.data,path = paste("./results/results_replicates_ID",ID,".csv",sep=""),delim = ";",append = F)
  # delete temporary files
  file.remove(paste("./results/",num.files,sep=""))
}
