
# palettes for empirical and simulated data
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

# in configurations 6,7, and 11, topology is always linear. useful to name the files
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
  
  data.list[[i.ID]]$net.sign <- as.factor(ifelse(data.list[[i.ID]]$net.effect>0,"positive","negative"))
  data.list[[i.ID]]$direct.sign <- as.factor(ifelse(data.list[[i.ID]]$direct.effect>0,"positive","negative"))
  
  data.list[[i.ID]]$patch <- ifelse(data.list[[i.ID]]$spatial.distance == 0,"within","across")
  
  # extract some information from the raw data, prior to averaging
  sim.ratios <- data.list[[i.ID]] %>% group_by(simulation.ID,patch,net.sign) %>% 
    summarise(number = n()) %>% spread(net.sign,number) %>% summarise(pos.neg.ratio = positive/negative)
  
  sim.magnitudes <- data.list[[i.ID]] %>% group_by(simulation.ID,patch,net.sign) %>% summarise(mean.effect = mean(net.effect),
                                                                               median.effect = median(net.effect),
                                                                               sd.effect = sd(net.effect)) 
  sim.switches <- data.list[[i.ID]] %>% filter(direct.effect != 0) %>% group_by(simulation.ID,patch,direct.sign) %>% summarise(count = n(),sign.switches = sum(sign(direct.effect) != sign(net.effect)))
  
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

###################################
###################################
# table 1

# 1 - +/- ratio: ratios
# readr::write_delim(x = ratios,path = paste("/home/david/CREAF/FPU/spatial_metaecosystem/results/ratios_table.csv",sep=""),delim = ";")
# # 2 - mean magnitude: magnitudes
# readr::write_delim(x = magnitudes,path = paste("/home/david/CREAF/FPU/spatial_metaecosystem/results/magnitudes_table.csv",sep=""),delim = ";")
# # 3 - sign switches: sign.switches
# readr::write_delim(x = sign.switches,path = paste("/home/david/CREAF/FPU/spatial_metaecosystem/results/switches_table.csv",sep=""),delim = ";")
# sign.switches %>% group_by(simulation.ID,patch,direct.sign) %>% summarize(freq.switch = sign.switches/count)

##########################################
##########################################
sp.effect.data <- full.data %>% group_by(replicate,trigger.sp,affected.sp,simulation.ID,patch) %>% summarise(mean.effect = mean(net.effect))

sp.effect.data <- spread(sp.effect.data,key = patch,value = mean.effect)
sp.effect.data <- droplevels(subset(sp.effect.data,affected.sp != 0))
sp.effect.data$affected.sp <- as.integer(sp.effect.data$affected.sp)
sp.effect.data$pairs <- as.factor(ifelse(sp.effect.data$affected.sp==sp.effect.data$trigger.sp,"i==j","i!=j"))

# in case I want to generate a contour around the set of points
contour.data <- NULL
prob <- c(0.8)
n.points <- 500#2000

for(i.sim in 1:length(unique(sp.effect.data$simulation.ID))){
  for(i.pairs in 1:length(unique(sp.effect.data$pairs))){
    my.data <- subset(sp.effect.data,
                      simulation.ID == unique(sp.effect.data$simulation.ID)[i.sim] &
                      pairs == unique(sp.effect.data$pairs)[i.pairs])

    kernel.data <- kde2d(my.data$within,my.data$across,n = n.points)
    dx <- diff(kernel.data$x[1:2])  # lifted from emdbook::HPDregionplot()
    dy <- diff(kernel.data$y[1:2])
    sz <- sort(kernel.data$z)
    c1 <- cumsum(sz) * dx * dy

    # plot:
    dimnames(kernel.data$z) <- list(kernel.data$x,kernel.data$y)
    dc <- reshape2::melt(kernel.data$z)
    dc$prob <- approx(sz,1-c1,dc$value)$y
    dc$simulation.ID <- unique(sp.effect.data$simulation.ID)[i.sim]
    dc$pairs <- unique(sp.effect.data$pairs)[i.pairs]
    contour.data <- rbind(contour.data,dc)
  }
}# for i.sim

# plot
cuadrants.plot <- ggplot(contour.data,aes(x = Var1,y = Var2))
cuadrants.plot <- cuadrants.plot +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "grey80",alpha = 0.7) +
  annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "grey80",alpha = 0.7) +
  geom_hline(yintercept = 0,color="grey60") + 
  geom_vline(xintercept = 0,color="grey60") + 
  # geom_abline(intercept=0, slope=1,color="grey70") +
  geom_point(aes(x=within,y=across,color=pairs),data=sp.effect.data,size = 0.8,alpha = 0.2) +
  # geom_contour(aes(z=prob,linetype=pairs),color = "black",size = 0.2,breaks=prob)+#,linetype = "dashed") +
  
  scale_color_manual(name = "",
                     labels = c("intraspecific effect","interspecific effect"),
                     values = cbPalette[c(6,7)]) +
  # scale_color_hue(l=40) +
  guides(colour = guide_legend(override.aes = list(alpha=1,size = 3)),linetype=F,fill=F) +
  facet_grid(simulation.ID~.) +
  theme_Publication() + 
  theme(panel.border = element_rect(colour="black")) +
  xlab("intra-patch net effect") + ylab("inter-patch net effect") +
  scale_x_continuous(breaks = seq(-1.6,1.2,0.2)) + ylim(-0.3,0.2) + xlim(-1.6,1.2)

tiff(paste("./results/images/combined_effects_",topology,".tiff",sep=""), res=600, compression = "lzw", width = 4500, height = 3500, units = "px")
print(cuadrants.plot)
dev.off()

##################################
##################################
# spatial distance and path length line plots
p_dodge <- position_dodge(0.2)

distances.data <- full.data %>% filter(path.length != 0) %>% group_by(simulation.ID,spatial.distance,path.length) %>% summarise(avg.net.effect = mean(abs(net.effect)),
                                                                                                                                                      lower.error = quantile(abs(net.effect),0.025),
                                                                                                                                                      upper.error = quantile(abs(net.effect),0.975))

distances.plot <- ggplot(distances.data,aes(x = spatial.distance,y = avg.net.effect, group = simulation.ID))

distances.plot <- distances.plot + stat_summary(aes(color = simulation.ID),
                                                    fun.y = mean,
                                                    # fun.ymin = function(x) ifelse(mean(x) - sd(x) > 0,mean(x) - sd(x),0), 
                                                    fun.ymin = function(x) mean(x) - quantile(x,0.025),
                                                    fun.ymax = function(x) mean(x) + quantile(x,0.975), 
                                                    geom = "errorbar",
                                                    position = p_dodge, width = 0.4) 

distances.plot <- distances.plot + stat_summary(aes(color = simulation.ID),
                                                    fun.y = mean, 
                                                    geom = "line", 
                                                    position = p_dodge, 
                                                    size = 1.5)

distances.plot <- distances.plot + stat_summary(fun.y = mean,
                                                    geom = "point", 
                                                    position = p_dodge, 
                                                    size = 2.1, 
                                                    shape = 21, 
                                                    fill = "white")

distances.plot <- distances.plot + scale_color_manual(name = "", values = cbPalette[c(2,4,6)])
# distances.plot <- distances.plot + scale_color_grey(name = "",start = 0.8,end = 0.2)
distances.plot <- distances.plot + DGC::theme_Publication() 
distances.plot <- distances.plot + xlab("spatial distance") + ylab("net effect") 
# distances.plot <- distances.plot + ylim(0,0.26)
distances.plot <- distances.plot + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1.1,0.1))

tiff(paste("./results/images/spatial_distance_",topology,".tiff",sep=""), res=600, compression = "lzw", width = 4000, height = 2500, units = "px")
print(distances.plot)
dev.off()

###############
###############
# averaged by path length

path.length.plot <- ggplot(distances.data,aes(x = path.length,y = avg.net.effect, group = simulation.ID))

path.length.plot <- path.length.plot + stat_summary(aes(color = simulation.ID),
                                                fun.y = mean,
                                                # fun.ymin = function(x) ifelse(mean(x) - sd(x) > 0,mean(x) - sd(x),0), 
                                                # fun.ymax = function(x) mean(x) + sd(x), 
                                                fun.ymin = function(x) mean(x) - quantile(x,0.025),
                                                fun.ymax = function(x) mean(x) + quantile(x,0.975), 
                                                geom = "errorbar",
                                                position = p_dodge, width = 0.4) 

path.length.plot <- path.length.plot + stat_summary(aes(color = simulation.ID),
                                                fun.y = mean, 
                                                geom = "line", 
                                                position = p_dodge, 
                                                size = 1.5)

path.length.plot <- path.length.plot + stat_summary(fun.y = mean,
                                                geom = "point", 
                                                position = p_dodge, 
                                                size = 2.1, 
                                                shape = 21, 
                                                fill = "white")

path.length.plot <- path.length.plot + scale_color_manual(name = "", values = cbPalette[c(2,4,6)])
# path.length.plot <- path.length.plot + scale_color_grey(name = "",start = 0.8,end = 0.2)
path.length.plot <- path.length.plot + DGC::theme_Publication()
path.length.plot <- path.length.plot + xlab("average path length") + ylab("net effect") 
path.length.plot <- path.length.plot + scale_x_continuous(limits = c(0.85,7.2),breaks = 1:7)
path.length.plot <- path.length.plot + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1.1,0.1))
# path.length.plot <- path.length.plot + ylim(0,0.7)

tiff(paste("./results/images/path_length_",topology,".tiff",sep=""), res=600, compression = "lzw", width = 4300, height = 2500, units = "px")
print(path.length.plot)
dev.off()

