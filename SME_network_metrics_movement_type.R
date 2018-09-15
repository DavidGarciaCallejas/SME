
# palettes for empirical and simulated data
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#########################
# generate_label_df <- function(my.data, HSD, flev, vspace = 0.5, my.var){
#   Tukey.levels <- HSD[[flev]][, 4]
#   Tukey.labels <- multcompView::multcompLetters(Tukey.levels)["Letters"]
#   plot.labels <- names(Tukey.labels[["Letters"]])
#   boxplot.df <- plyr::ddply(my.data, flev, function(x) max(fivenum(my.data[,my.var])) + 
#                               vspace)
#   plot.levels <- data.frame(plot.labels, labels = Tukey.labels[["Letters"]], 
#                             stringsAsFactors = FALSE)
#   labels.df <- merge(plot.levels, boxplot.df, by.x = "plot.labels", 
#                      by.y = flev, sort = FALSE)
#   return(labels.df)
# }
##########################
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
  
  # data.list[[i.ID]] <- readr::read_delim(paste("./results/results_replicates_ID",ID[i.ID],".csv",sep=""),delim = ";",col_types = cols())
  # temp.species.data <- readr::read_delim(paste("./results/species_data_replicates_ID",ID[i.ID],".csv",sep=""),delim = ";",col_types = cols())
  
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
# full.data$sign <- as.factor(ifelse(full.data$net.effect>0,"positive","negative"))
# full.data$sign <- as.factor(ifelse(full.data$avg.net.effect>0,"positive","negative"))

# full.data$patch <- ifelse(full.data$spatial.distance == 0,"within","across")

rm(data.list)
gc()

######################
######################
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

# 1 - connectance
connectance.results <- expand.grid(replicate = unique(full.data$replicate),
                                   simulation = simulation.levels,
                                   patch = patch)
connectance.results$connectance <- 0

for(i.sim in 1:length(simulation.levels)){
  
  sim.data <- subset(full.data,simulation.ID == simulation.levels[i.sim])
  
  for(i.rep in 1:replicates){
    rep.data <- subset(sim.data,replicate == i.rep)
    interpatch.data <- subset(rep.data,patch == "across")
    intrapatch.data <- subset(rep.data,patch == "within")
    num.communities <- rep.data$N[1]
    num.sp <- rep.data$S[1]
    
    # intrapatch connectance
    # denominator is multiplied by N because I have all the interactions of all patches together, 
    # and I did not include a column with the specific patch in which the interaction takes place
    connectance.results$connectance[connectance.results$replicate == i.rep & 
                                      connectance.results$simulation == simulation.levels[i.sim] & 
                                      connectance.results$patch == "intra-patch"] <- sum(intrapatch.data$direct.effect != 0)/((num.sp^2)*num.communities)
    # interpatch connectance
    # the set of (directed) potential links is S*(S*(N-1))*N
    potential.interpatch.links <- num.sp * (num.sp * (num.communities - 1)) * num.communities
    connectance.results$connectance[connectance.results$replicate == i.rep & 
                                      connectance.results$simulation == simulation.levels[i.sim] & 
                                      connectance.results$patch == "inter-patch"] <- sum(interpatch.data$direct.effect != 0)/potential.interpatch.links
  }# for i.rep
}# for i.sim

connectance.spread <- spread(connectance.results,key = patch,value = connectance)
names(connectance.spread)[3] <- "intra-patch.connectance"
names(connectance.spread)[4] <- "inter-patch.connectance"

# kruskal wallis for differences in intra and interpatch connectances among groups

connectance.plot = ggplot(connectance.results) + 
  geom_boxplot(aes(x = simulation, y = connectance, fill = patch))

################
# other metrics

network.metrics <- expand.grid(replicate = unique(full.data$replicate),
                              simulation = simulation.levels)
network.metrics <- left_join(network.metrics,connectance.spread)
network.metrics$in.degree <- 0
network.metrics$out.degree <- 0
network.metrics$clustering.coefficient <- 0
network.metrics$average.path.length <- 0
network.metrics$quantitative.link.density <- 0
network.metrics$quantitative.connectance <- 0

for(i.sim in 1:length(simulation.levels)){
  
  sim.data <- subset(direct.edges.data,simulation.ID == simulation.levels[i.sim])
  
  for(i.rep in 1:replicates){
    
    result.index <- which(network.metrics$replicate == i.rep & network.metrics$simulation == simulation.levels[i.sim])
    
    rep.data <- subset(sim.data,replicate == i.rep)
    rep.data <- rep.data[,c("trigger.sp.ID","affected.sp.ID","direct.effect")]
    
    # build igraph object and adjacency matrix
    my.network <- graph.data.frame(rep.data)
    E(my.network)$weight = rep.data$direct.effect
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
    
  }# for i.rep
}# for i.sim

network.metrics.gather <- gather(network.metrics,key = "metric",value = "value",-replicate,-simulation)
my.labels <- c("dispersal","foraging","dispersal\nand foraging")

##################
# apply here any subset for the variables to be shown
network.metrics.gather <- droplevels(subset(network.metrics.gather, metric %in% c("average.path.length","clustering.coefficient","quantitative.connectance")))

###################
# tukey test for statistical significance of differences between groups

path.length.test <- aov(average.path.length ~ simulation, data = network.metrics)
path.length.tukey <- TukeyHSD(path.length.test,ordered = FALSE, conf.level = 0.95)
path.length.labels <- generate_label_df(my.data = network.metrics,
                                        HSD = path.length.tukey,
                                        flev = "simulation",
                                        my.var = "average.path.length")
path.length.labels$metric <- "average.path.length"
#### hack - groups are all different but letters are not ordered, so do it manually
path.length.labels$labels <- c("b","c","a")

clustering.coefficient.test <- aov(clustering.coefficient ~ simulation, data = network.metrics)
clustering.coefficient.tukey <- TukeyHSD(clustering.coefficient.test,ordered = FALSE, conf.level = 0.95)
clustering.coefficient.labels <- generate_label_df(my.data = network.metrics,
                                                   HSD = clustering.coefficient.tukey,
                                                   flev = "simulation",
                                                   my.var = "clustering.coefficient")
clustering.coefficient.labels$metric <- "clustering.coefficient"
#### hack - groups are all different but letters are not ordered, so do it manually
clustering.coefficient.labels$labels <- c("b","c","a")

quant.connectance.test <- aov(quantitative.connectance ~ simulation, data = network.metrics)
quant.connectance.tukey <- TukeyHSD(quant.connectance.test,ordered = FALSE, conf.level = 0.95)
quant.connectance.labels <- generate_label_df(my.data = network.metrics,
                                        HSD = quant.connectance.tukey,
                                        flev = "simulation",
                                        my.var = "quantitative.connectance")
quant.connectance.labels$metric <- "quantitative.connectance"
#### hack - groups are all different but letters are not ordered, so do it manually
quant.connectance.labels$labels <- c("b","c","a")

metric.labels <- bind_rows(path.length.labels,
                           clustering.coefficient.labels,
                           quant.connectance.labels)

# again, modify manually the heights of the letters
metric.labels$V1[metric.labels$metric == "clustering.coefficient"] <- 0.7
metric.labels$V1[metric.labels$metric == "quantitative.connectance"] <- 0.1

##################
# final touches
labels <- c(average.path.length = "average path length", 
            clustering.coefficient = "clustering coefficient", 
            quantitative.connectance = "quantitative connectance")

##################
# plot
metrics.plot <- ggplot(network.metrics.gather) + 
  geom_boxplot(aes(x = simulation,y = value,fill = simulation)) +
  geom_text(data = metric.labels, aes(x = plot.labels, y = V1, label = labels)) +
  scale_fill_manual(values = cbPalette[c(2,4,6)]) + 
  scale_x_discrete(labels = my.labels) +
  xlab("") + ylab("metric value") +
  guides(fill = FALSE) +
  theme_Publication() +
  facet_wrap(~metric,scales = "free_y",labeller = labeller(metric = labels)) +
  NULL

readr::write_delim(network.metrics,path = paste("./results/",topology,"/network_metrics.csv",sep=""),delim = ";")
tiff(paste("./results/",topology,"/images/network_metrics_",topology,".tiff",sep=""), res=600, compression = "lzw", width = 4000, height = 2500, units = "px")
print(metrics.plot)
dev.off()
