####
# analysis of network structural metrics
# in this script, differentiating by movement mode and interaction type,
# as in Fig. 3 of the manuscript.
# In this case, no statistical tests are performed.
####

# palettes for empirical and simulated data
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##########################
# auxiliary functions for some metrics

hill.diversity <- function(abundances,q = 1){
  abundances <- abundances[abundances != 0]
  abundances <- abundances/sum(abundances)
  R <- length(abundances)
  # hill diversity is not defined for q = 1,
  # but its limit exists and equals
  # the exponential of shannon entropy 
  # (Jost 2006,2007,Tuomisto 2012)
  if(q == 1){
    D <- exp(-sum(abundances*log(abundances)))
  }else{
    mean.p <- (sum(abundances*abundances^(q-1)))^(1/(q-1))
    D <- 1/mean.p
  }
  return(D)
}

QuantitativeLinkDensity <- function(interaction.matrix){
  # after LDq' in Banasek-Richter et al. 2009
  
  # check consistency of the coefficients, they need to be > 0
  # for calculating the shannon index
  if(sum(interaction.matrix < 0) > 0){
    interaction.matrix <- abs(interaction.matrix)
  }
  
  # 1 - in and out degrees given by the effective number of links
  num.sp <- nrow(interaction.matrix)
  in.degrees <- apply(interaction.matrix,MARGIN = 2,FUN = hill.diversity)
  out.degrees <- apply(interaction.matrix,MARGIN = 1,FUN = hill.diversity)
  
  # 2 - link density
  LD <- 1/(2*num.sp) * (sum(in.degrees) + sum(out.degrees))
  
  LD
}

QuantitativeConnectance <- function(interaction.matrix){
  # after LDq' in Banasek-Richter et al. 2009
  
  # check consistency of the coefficients, they need to be > 0
  # for calculating the shannon index
  if(sum(interaction.matrix < 0) > 0){
    interaction.matrix <- abs(interaction.matrix)
  }
  
  # 1 - in and out degrees given by the effective number of links
  num.sp <- nrow(interaction.matrix)
  in.degrees <- apply(interaction.matrix,MARGIN = 2,FUN = hill.diversity)
  out.degrees <- apply(interaction.matrix,MARGIN = 1,FUN = hill.diversity)
  
  # 2 - link density
  LD <- 1/(2*num.sp) * (sum(in.degrees) + sum(out.degrees))
  
  # 3 - connectance
  (LD/num.sp)
}

#########################
#########################
# set of simulations
# 1 - global, constant dispersal
# 2 - linear, constant dispersal
# 3 - random, constant dispersal
# 4 - linear, constant foraging
# 5 - linear, scaling foraging
# 6 - as 2, with only antagonism, c = 0.2
# 7 - as 4, with only antagonism, c = 0.2
# 8 - as 5, with only antagonism, c = 0.2
#######################
# 9 - linear, constant dispersal,constant foraging
# 10 - linear, scaling dispersal, scaling foraging
# 11 - as 9, with only antagonism, c = 0.2
# 12 - as 10, with only antagonism, c = 0.2
ID <- c(6,7,11)
file.version <- 4

topology <- "linear"

##################
##################
# dispersal and foraging rates
# Other rates will be discarded
dispersal.rate <- 0.5
foraging.rate <- 0.5

pars = read.table("./data/parameters_set.csv",header = T,sep = ";",dec = ".",stringsAsFactors = F)
species.data <- NULL
data.list <- list()
ratios <- NULL
magnitudes <- NULL
sign.switches <- NULL

for(i.ID in 1:length(ID)){
  
  data.list[[i.ID]] <- readr::read_delim(paste("./results/results_replicates_ID",ID[i.ID],".csv",sep=""),delim = ";",col_types = cols())
  temp.species.data <- readr::read_delim(paste("./results/species_data_replicates_ID",ID[i.ID],".csv",sep=""),delim = ";",col_types = cols())
  
  # trophic level of the affected species
  data.list[[i.ID]] <- left_join(data.list[[i.ID]],temp.species.data,by=c("replicate", "affected.sp" = "species"))
  names(data.list[[i.ID]])[names(data.list[[i.ID]]) %in% c("PreyAvg.TL","niche.axis")] <- c("affected.sp.TL","affected.sp.niche")
  
  # trophic level of the source species
  data.list[[i.ID]] <- left_join(data.list[[i.ID]],temp.species.data,by=c("replicate", "trigger.sp" = "species"))
  names(data.list[[i.ID]])[names(data.list[[i.ID]]) %in% c("PreyAvg.TL","niche.axis")] <- c("trigger.sp.TL","trigger.sp.niche")
  
  data.list[[i.ID]]$sign <- as.factor(ifelse(data.list[[i.ID]]$net.effect>0,"positive","negative"))
  data.list[[i.ID]]$patch <- ifelse(data.list[[i.ID]]$spatial.distance == 0,"within","across")
  
  # extract some information from the raw data, prior to averaging
  sim.ratios <- data.list[[i.ID]] %>% group_by(simulation.ID,patch,sign) %>% 
    summarise(number = n()) %>% spread(sign,number) %>% summarise(pos.neg.ratio = positive/negative)
  
  sim.magnitudes <- data.list[[i.ID]] %>% group_by(simulation.ID,patch,sign) %>% summarise(mean.effect = mean(net.effect),
                                                                                           median.effect = median(net.effect),
                                                                                           sd.effect = sd(net.effect)) 
  sim.switches <- data.list[[i.ID]] %>% filter(direct.effect != 0) %>% group_by(simulation.ID,patch,sign) %>% summarise(count = n(),sign.switches = sum(sign(direct.effect) != sign(net.effect)))
  
  ratios <- rbind(ratios,sim.ratios)
  magnitudes <- rbind(magnitudes,sim.magnitudes)
  sign.switches <- rbind(sign.switches,sim.switches)
  
  data.list[[i.ID]]$topology <- pars$topology[ID[i.ID]]
  data.list[[i.ID]]$N <- pars$N[ID[i.ID]]
  data.list[[i.ID]]$S <- pars$S[ID[i.ID]]
  data.list[[i.ID]] <- data.list[[i.ID]][complete.cases(data.list[[i.ID]]),]
  temp.species.data$simulation.ID <- ID[i.ID]
  species.data <- rbind(species.data,temp.species.data)
}

full.data <- dplyr::bind_rows(data.list)
full.data$simulation.ID <- as.factor(full.data$simulation.ID)
full.data$simulation.ID <- plyr::revalue(full.data$simulation.ID,replace = c("1" = "global-const disp",
                                                                             "2" = "linear-const disp",
                                                                             "3" = "random-const disp",
                                                                             "4" = "linear-const forag",
                                                                             "5" = "linear-scaling forag",
                                                                             "6" = "dispersal",
                                                                             "7" = "foraging",
                                                                             "8" = "linear-scaling forag-antag",
                                                                             "9" = "linear-const disp-const forag",
                                                                             "10" = "linear-scaling disp-scaling forag",
                                                                             "11" = "dispersal and foraging",
                                                                             "12" = "linear-scaling disp-scaling forag-antag"))
species.data$simulation.ID <- as.factor(species.data$simulation.ID)
species.data$simulation.ID <- plyr::revalue(species.data$simulation.ID,replace = c("1" = "global-const disp",
                                                                                   "2" = "linear-const disp",
                                                                                   "3" = "random-const disp",
                                                                                   "4" = "linear-const forag",
                                                                                   "5" = "linear-scaling forag",
                                                                                   "6" = "dispersal",
                                                                                   "7" = "foraging",
                                                                                   "8" = "linear-scaling forag-antag",
                                                                                   "9" = "linear-const disp-const forag",
                                                                                   "10" = "linear-scaling disp-scaling forag",
                                                                                   "11" = "dispersal and foraging",
                                                                                   "12" = "linear-scaling disp-scaling forag-antag"))

rm(data.list)
gc()

######################
######################
effect.levels <- c("direct","net")
simulation.levels <- unique(full.data$simulation.ID)
patch <- c("intra-patch","inter-patch")
replicates <- max(full.data$replicate)

######################
######################
# network metrics

# 0 - data preparation
full.data$trigger.sp.ID <- paste(full.data$trigger.sp,full.data$trigger.sp.patch,sep="_")
full.data$affected.sp.ID <- paste(full.data$affected.sp,full.data$affected.sp.patch,sep="_")

direct.edges.data <- subset(full.data,direct.effect != 0)
direct.edges.data <- direct.edges.data[,c("simulation.ID","replicate","trigger.sp.ID","affected.sp.ID","direct.effect")]

net.edges.data <- subset(full.data,net.effect != 0)
net.edges.data <- net.edges.data[,c("simulation.ID","replicate","trigger.sp.ID","affected.sp.ID","net.effect")]

# 1 - connectance
connectance.results <- expand.grid(replicate = unique(full.data$replicate),
                                   effect = effect.levels,
                                   simulation = simulation.levels,
                                   patch = patch)
connectance.results$connectance <- 0

for(i.effect in 1:length(effect.levels)){
  for(i.sim in 1:length(simulation.levels)){
    
    sim.data <- subset(full.data,simulation.ID == simulation.levels[i.sim])
    
    for(i.rep in 1:replicates){
      rep.data <- subset(sim.data,replicate == i.rep)
      interpatch.data <- subset(rep.data,patch == "across")
      intrapatch.data <- subset(rep.data,patch == "within")
      num.communities <- rep.data$N[1]
      num.sp <- rep.data$S[1]
      
      my.intra.num <- 0
      my.inter.num <- 0
      if(effect.levels[i.effect] == "direct"){
        my.intra.num <- sum(intrapatch.data$direct.effect != 0)
        my.inter.num <- sum(interpatch.data$direct.effect != 0)
      }else{
        my.intra.num <- sum(intrapatch.data$net.effect != 0)
        my.inter.num <- sum(interpatch.data$net.effect != 0)
      }
      
      # intrapatch connectance
      # denominator is multiplied by N because I have all the interactions of all patches together, 
      # and I did not include a column with the specific patch in which the interaction takes place
      connectance.results$connectance[connectance.results$replicate == i.rep & 
                                        connectance.results$effect == effect.levels[i.effect] &
                                        connectance.results$simulation == simulation.levels[i.sim] & 
                                        connectance.results$patch == "intra-patch"] <- my.intra.num/((num.sp^2)*num.communities)
      # interpatch connectance
      # the set of (directed) potential links is S*(S*(N-1))*N
      potential.interpatch.links <- num.sp * (num.sp * (num.communities - 1)) * num.communities
      connectance.results$connectance[connectance.results$replicate == i.rep & 
                                        connectance.results$effect == effect.levels[i.effect] &
                                        connectance.results$simulation == simulation.levels[i.sim] & 
                                        connectance.results$patch == "inter-patch"] <- my.inter.num/potential.interpatch.links
    }# for i.rep
  }# for i.sim
}# for i.effect

connectance.spread <- spread(connectance.results,key = patch,value = connectance)
names(connectance.spread)[4] <- "intra-patch.connectance"
names(connectance.spread)[5] <- "inter-patch.connectance"

# take a look at the connectances
# connectance.plot = ggplot(connectance.results) + 
#   geom_boxplot(aes(x = simulation, y = connectance, fill = patch)) + 
#   facet_grid(.~effect) +
#   NULL
# connectance.plot

################
# other metrics

network.metrics <- expand.grid(replicate = unique(full.data$replicate),
                               effect = effect.levels,
                              simulation = simulation.levels)

network.metrics <- left_join(network.metrics,connectance.spread)
network.metrics$in.degree <- 0
network.metrics$out.degree <- 0
network.metrics$clustering.coefficient <- 0
network.metrics$average.path.length <- 0
network.metrics$quantitative.link.density <- 0
network.metrics$quantitative.connectance <- 0
network.metrics$quantitative.modularity <- 0
network.metrics$number.modules <- 0

# 
for(i.sim in 1:length(simulation.levels)){
  for(i.effect in 1:length(effect.levels)){
    
    if(effect.levels[i.effect] == "direct"){
      sim.data <- subset(direct.edges.data,simulation.ID == simulation.levels[i.sim])
    }else{
      sim.data <- subset(net.edges.data,simulation.ID == simulation.levels[i.sim])
    }
    
    for(i.rep in 1:replicates){
      
      result.index <- which(network.metrics$replicate == i.rep & 
                              network.metrics$simulation == simulation.levels[i.sim] &
                              network.metrics$effect == effect.levels[i.effect])
      
      rep.data <- subset(sim.data,replicate == i.rep)
      rep.data <- rep.data[,c(3,4,5)]
      names(rep.data)[3] <- "effect"
      
      # build igraph object and adjacency matrix
      my.network <- graph.data.frame(rep.data)
      E(my.network)$weight = rep.data$effect
      my.interaction.matrix <- as_adjacency_matrix(my.network,attr = "weight")
      
      # degree
      network.metrics$in.degree[result.index] <- mean(degree(my.network,mode = "in"))
      network.metrics$out.degree[result.index] <- mean(degree(my.network,mode = "out"))
      
      # clustering coefficient
      network.metrics$clustering.coefficient[result.index] <- transitivity(my.network,type = "global")
      
      # average path length
      network.metrics$average.path.length[result.index] <- mean_distance(my.network,directed = TRUE)
      
      # quantitative link density
      network.metrics$quantitative.link.density[result.index] <- QuantitativeLinkDensity(my.interaction.matrix)
      
      # quantitative connectance
      network.metrics$quantitative.connectance[result.index] <- QuantitativeConnectance(my.interaction.matrix)
      
      # quantitative modularity -- including negative weights!! see the function help page
      # preliminary tests showed that no more than 10 modules are usually found, so take that number as a maximum
      # in order to make it a bit faster
      my.community <- igraph::cluster_spinglass(my.network,implementation = "neg", spins = 10)
      
      network.metrics$quantitative.modularity[result.index] <- my.community$modularity
      network.metrics$number.modules[result.index] <- length(unique(my.community$membership))
      
    }# for i.rep
  }# for i.effect
}# for i.sim
# 

# prepare the data for the image
network.metrics.gather <- gather(network.metrics,key = "metric",value = "value",-replicate,-simulation,-effect)

# also, take a look at mean values
mean.metrics <- network.metrics.gather %>% group_by(effect,simulation,metric) %>% summarise(mean.value = mean(value), sd.value = sd(value))
mean.metrics <- subset(mean.metrics, metric %in% c("inter-patch.connectance","intra-patch.connectance","average.path.length","quantitative.modularity"))

##################
# apply here any subset for the variables to be shown
network.metrics.plot <- droplevels(subset(network.metrics.gather, metric %in% c("average.path.length",
                                                                                "intra-patch.connectance",
                                                                                "inter-patch.connectance",
                                                                                  # "clustering.coefficient",
                                                                                  # "quantitative.connectance",
                                                                                  "quantitative.modularity")))

# order
network.metrics.plot$metric <- factor(network.metrics.plot$metric, levels = c("intra-patch.connectance",
                                                                                 "inter-patch.connectance",
                                                                                 "average.path.length",
                                                                                 "quantitative.modularity"))
# and write nice names
network.metrics.plot$metric <- plyr::revalue(network.metrics.plot$metric, c("average.path.length" = "average path length",
                                                                            "intra-patch.connectance" = "intra-patch connectance",
                                                                            "inter-patch.connectance" = "inter-patch connectance",
                                                                            # clustering.coefficient = "clustering coefficient",
                                                                            # quantitative.connectance = "quantitative connectance",
                                                                            "quantitative.modularity" = "modularity"))

sim.labels <- c("dispersal","foraging","dispersal\nand foraging")

##################
# plot

metrics.plot <- 
  ggplot(network.metrics.plot, aes(simulation, value, fill = effect)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge", color = "grey20") +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, 
               width = 0.3,
               position = position_dodge(.9),
               color = "grey20")+
  facet_wrap(~metric,scales = "free_y") +
  # scale_fill_manual(values = c("grey30","grey70"))+
  # scale_fill_manual(values = cbPalette[c(5,7)])+
  scale_fill_manual(values = c("darkolivegreen","darkolivegreen3"))+
  scale_x_discrete(labels = sim.labels) +
  xlab("") + ylab("metric value") +
  DGC::theme_Publication()+
  theme(strip.background = element_blank()) +
  guides(fill = FALSE)+
  NULL

# other types of plot are not that nice, 
# considering that e.g. connectance is always 1 for some networks

# metrics.plot <- ggplot(network.metrics.gather) + 
  # geom_boxplot(aes(x = simulation,y = value,fill = effect)) +
  # scale_fill_manual(values = cbPalette[c(2,4,6)]) + 
  # scale_x_discrete(labels = sim.labels) +
  # xlab("") + ylab("metric value") +
  # guides(fill = FALSE) +
  # DGC::theme_Publication() +
  # facet_wrap(~metric,scales = "free_y")+#,labeller = labeller(metric = metric.labels)) +
  # NULL

##########
# store the figure and the metric results
# readr::write_delim(network.metrics,path = paste("./results/network_metrics.csv",sep=""),delim = ";")
tiff(paste("./results/images/network_metrics_",topology,".tiff",sep=""), res=600, compression = "lzw", width = 4000, height = 4000, units = "px")
print(metrics.plot)
dev.off()
