###################
# plot results of the sensitivity analyses
###################

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

ID <- c("6","6_low_disp_no_for","6_high_disp_no_for",
        "7","7_no_disp_low_for","7_no_disp_high_for",
        "11","11_low_disp_low_for","11_low_disp_high_for","11_high_disp_low_for","11_high_disp_high_for")
# ID <- ID[1:5]
file.version <- 5
topology <- "linear"

##################
##################
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
  
  id.num <- as.numeric(substr(ID[i.ID],1,1))

  data.list[[i.ID]]$simulation.ID <- as.character(data.list[[i.ID]]$simulation.ID)
  data.list[[i.ID]]$topology <- pars$topology[id.num]
  data.list[[i.ID]]$N <- pars$N[id.num]
  data.list[[i.ID]]$S <- pars$S[id.num]
  data.list[[i.ID]] <- data.list[[i.ID]][complete.cases(data.list[[i.ID]]),]
  temp.species.data$simulation.ID <- ID[i.ID]
  species.data <- rbind(species.data,temp.species.data)
}

full.data <- dplyr::bind_rows(data.list)
full.data$simulation.ID <- as.factor(full.data$simulation.ID)
full.data$simulation.ID <- plyr::revalue(full.data$simulation.ID,replace = c("6_low_disp_no_for" = "dispersal - low d",
                                                                             "6" = "dispersal - intermediate d",
                                                                             "6_high_disp_no_for" = "dispersal - high d",
                                                                             "7_no_disp_low_for" = "foraging - low f",
                                                                             "7" = "foraging - intermediate f",
                                                                             "7_no_disp_high_for" = "foraging - high f",
                                                                             "11" = "combined - intermediate d, intermediate f",
                                                                             "11_low_disp_low_for" = "combined - low d, low f",
                                                                             "11_low_disp_high_for" = "combined - low d, high f",
                                                                             "11_high_disp_low_for" = "combined - high d, low f",
                                                                             "11_high_disp_high_for" = "combined - high d, high f"))
species.data$simulation.ID <- as.factor(species.data$simulation.ID)
species.data$simulation.ID <- plyr::revalue(species.data$simulation.ID,replace = c("6_low_disp_no_for" = "dispersal - low d",
                                                                                   "6" = "dispersal - intermediate d",
                                                                                   "6_high_disp_no_for" = "dispersal - high d",
                                                                                   "7_no_disp_low_for" = "foraging - low f",
                                                                                   "7" = "foraging - intermediate f",
                                                                                   "7_no_disp_high_for" = "foraging - high f",
                                                                                   "11" = "combined - intermediate d, intermediate f",
                                                                                   "11_low_disp_low_for" = "combined - low d, low f",
                                                                                   "11_low_disp_high_for" = "combined - low d, high f",
                                                                                   "11_high_disp_low_for" = "combined - high d, low f",
                                                                                   "11_high_disp_high_for" = "combined - high d, high f"))

my.sims <- unique(full.data$simulation.ID)
full.data$simulation.ID <- as.character(full.data$simulation.ID)
full.data$dispersal.param <- "0"
full.data$foraging.param <- "0"

for(i.sim in 1:length(unique(full.data$simulation.ID))){
  my.index <- which(full.data$simulation.ID == my.sims[i.sim])
  if(grepl(pattern = "low d",x = my.sims[i.sim])){
    full.data$dispersal.param[my.index] <- "low"
  }else if(grepl(pattern = "intermediate d",x = my.sims[i.sim])){
    full.data$dispersal.param[my.index] <- "intermediate"
  }else if(grepl(pattern = "high d",x = my.sims[i.sim])){
    full.data$dispersal.param[my.index] <- "high"
  }
  
  if(grepl(pattern = "low f",x = my.sims[i.sim])){
    full.data$foraging.param[my.index] <- "low"
  }else if(grepl(pattern = "intermediate f",x = my.sims[i.sim])){
    full.data$foraging.param[my.index] <- "intermediate"
  }else if(grepl(pattern = "high f",x = my.sims[i.sim])){
    full.data$foraging.param[my.index] <- "high"
  }
  
  if(grepl(pattern = "dispersal",x = my.sims[i.sim])){
    full.data$simulation.ID[my.index] <- "dispersal"
  }else if(grepl(pattern = "foraging",x = my.sims[i.sim])){
    full.data$simulation.ID[my.index] <- "foraging"
  }else if(grepl(pattern = "combined",x = my.sims[i.sim])){
    full.data$simulation.ID[my.index] <- "dispersal and foraging"
  }
}

full.data$patch[full.data$patch == "within"] <- "intra-patch"
full.data$patch[full.data$patch == "across"] <- "inter-patch"

rm(data.list)
gc()

###################################
###################################
# table 1

# extract some information from the raw data, prior to averaging
ratios <- full.data %>% group_by(simulation.ID,dispersal.param,foraging.param,patch,net.sign) %>% 
  summarise(number = n()) %>% spread(net.sign,number) %>% summarise(pos.neg.ratio = positive/negative)

magnitudes <- full.data %>% group_by(simulation.ID,dispersal.param,foraging.param,patch,net.sign) %>% summarise(mean.effect = mean(net.effect),
                                                                                         median.effect = median(net.effect),
                                                                                         sd.effect = sd(net.effect)) 
switches <- full.data %>% filter(direct.effect != 0) %>% group_by(simulation.ID,dispersal.param,foraging.param,patch,direct.sign) %>% summarise(count = n(),sign.switches = sum(sign(direct.effect) != sign(net.effect))) 

switches <- switches %>% group_by(simulation.ID,dispersal.param,foraging.param,patch,direct.sign) %>% summarise(switch.ratio = sign.switches/count)

###################################
###################################
# plot variations in table 1 across simulations

# tidy the data

my.sims <- unique(full.data$simulation.ID)

ratios$simulation.ID <- factor(ratios$simulation.ID, levels = c("dispersal", "foraging", "dispersal and foraging"))
levels(ratios$simulation.ID) <- gsub(" ", "\n", levels(ratios$simulation.ID))

magnitudes$simulation.ID <- factor(magnitudes$simulation.ID, levels = c("dispersal", "foraging", "dispersal and foraging"))
levels(magnitudes$simulation.ID) <- gsub(" ", "\n", levels(magnitudes$simulation.ID))

switches$simulation.ID <- factor(switches$simulation.ID, levels = c("dispersal", "foraging", "dispersal and foraging"))
levels(switches$simulation.ID) <- gsub(" ", "\n", levels(switches$simulation.ID))

magnitude.means <- magnitudes %>% filter(dispersal.param == "intermediate" | foraging.param == "intermediate") %>% 
  group_by(simulation.ID, patch) %>% summarise(mean.mag = mean(mean.effect))

switches.means <- switches %>% filter(dispersal.param == "intermediate" | foraging.param == "intermediate") %>% 
  group_by(simulation.ID, patch) %>% summarise(mean.sw = mean(switch.ratio))

#################
#################
# plot

test.ratios <- ggplot(ratios,aes(x = simulation.ID, group = interaction(simulation.ID,patch), fill = patch)) + 
  geom_boxplot(aes(y = pos.neg.ratio),outlier.alpha = 0.1) + 
  geom_point(data = ratios[ratios$dispersal.param == "intermediate" | ratios$foraging.param == "intermediate",],
             aes(x = simulation.ID, y = pos.neg.ratio, group = patch), fill = "darkgrey",size = 2, shape = 21,position = position_dodge(.75)) +
  scale_fill_manual(name = "",values = cbPalette[c(2,4)]) + 
  # theme_Publication() + 
  xlab("") + ylab("+/- ratio")
# test.ratios

test.magnitudes <- ggplot(magnitudes,aes(x = simulation.ID, group = interaction(simulation.ID,patch), fill = patch)) + 
  geom_boxplot(aes(y = mean.effect),outlier.size = 0.5) + 
  geom_point(data = magnitude.means,
             aes(x = simulation.ID, y = mean.mag, group = patch), fill = "darkgrey", size = 2, shape = 21,position = position_dodge(.75)) +
  scale_fill_manual(name = "",values = cbPalette[c(2,4)]) + 
  # theme_Publication() + 
  xlab("") + ylab("effect magnitude") + guides(fill = F)
# test.magnitudes

test.switches <- ggplot(switches,aes(x = simulation.ID, group = interaction(simulation.ID,patch), fill = patch)) + 
  geom_boxplot(aes(y = switch.ratio),outlier.size = 0.5) + 
  geom_point(data = switches.means,
             aes(x = simulation.ID, y = mean.sw, group = patch), fill = "darkgrey", size = 2, shape = 21,position = position_dodge(.75)) +
  scale_fill_manual(name = "",values = cbPalette[c(2,4)]) + 
  # theme_Publication() + 
  xlab("") + ylab("sign switch frequency") + guides(fill = F)
# test.switches

tiff(paste("./results/images/Fig_S1_v",file.version,".tiff",sep=""), res=600, compression = "lzw", width = 7000, height = 2500, units = "px")
print(test.ratios + test.magnitudes + test.switches)
dev.off()

##########################################
##########################################

full.data$dispersal.param <- paste(full.data$dispersal.param," d",sep="")
full.data$foraging.param <- paste(full.data$foraging.param," f",sep="")
full.data$parameter.values <- paste(full.data$dispersal.param,full.data$foraging.param,sep="-")
full.data$simulation.ID <- factor(full.data$simulation.ID,levels = c("dispersal","foraging","dispersal and foraging"))

full.data$parameter.values <- plyr::revalue(full.data$parameter.values,replace = c("low d-0 f" = "low d",
                                                                             "intermediate d-0 f" = "intermediate d",
                                                                             "high d-0 f" = "high d",
                                                                             "0 d-low f" = "low f",
                                                                             "0 d-intermediate f" = "intermediate f",
                                                                             "0 d-high f" = "high f"))
                                                                             # "11" = "combined - intermediate d, intermediate f",
                                                                             # "11_low_disp_low_for" = "combined - low d, low f",
                                                                             # "11_low_disp_high_for" = "combined - low d, high f",
                                                                             # "11_high_disp_low_for" = "combined - high d, low f",
                                                                             # "11_high_disp_high_for" = "combined - high d, high f"))

full.data$parameter.values <- factor(full.data$parameter.values,levels = c("low d-low f",
                                                                           "intermediate d-intermediate f", 
                                                                           "high d-high f", 
                                                                           "high d-low f", 
                                                                           "low d-high f",
                                                                           "low d","intermediate d", "high d",
                                                                           "low f","intermediate f", "high f"))

sp.effect.data <- full.data %>% group_by(replicate,trigger.sp,affected.sp,simulation.ID,parameter.values,patch) %>% summarise(mean.effect = mean(net.effect))

sp.effect.data <- spread(sp.effect.data,key = patch,value = mean.effect)
sp.effect.data <- droplevels(subset(sp.effect.data,affected.sp != 0))
sp.effect.data$affected.sp <- as.integer(sp.effect.data$affected.sp)
sp.effect.data$pairs <- as.factor(ifelse(sp.effect.data$affected.sp==sp.effect.data$trigger.sp,"i==j","i!=j"))

# cuadrants.plot <- ggplot(contour.data,aes(x = Var1,y = Var2))
cuadrants.plot <- ggplot()
cuadrants.plot <- cuadrants.plot +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "grey80",alpha = 0.7) +
  annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "grey80",alpha = 0.7) +
  geom_hline(yintercept = 0,color="grey60") + 
  geom_vline(xintercept = 0,color="grey60") + 
  # geom_abline(intercept=0, slope=1,color="grey70") +
  geom_point(aes(x=`intra-patch`,y=`inter-patch`,color=pairs),data=sp.effect.data,size = 0.8,alpha = 0.2) +
  # geom_contour(aes(z=prob,linetype=pairs),color = "black",size = 0.2,breaks=prob)+#,linetype = "dashed") +
  
  scale_color_manual(name = "",
                     labels = c("intraspecific effect","interspecific effect"),
                     values = cbPalette[c(6,7)]) +
  # scale_color_hue(l=40) +
  guides(colour = guide_legend(override.aes = list(alpha=1,size = 3)),linetype=F,fill=F) +
  facet_wrap(simulation.ID~parameter.values,drop = TRUE,ncol = 3) +
  # theme_Publication() + 
  theme(panel.border = element_rect(colour="black")) +
  xlab("intra-patch net effect") + ylab("inter-patch net effect") +
  scale_x_continuous(breaks = seq(-1.6,1.2,0.2)) + ylim(-0.3,0.2) + xlim(-1.6,1.2)

tiff(paste("./results/images/Fig_S2_v",file.version,".tiff",sep=""), res=600, compression = "lzw", width = 6000, height = 5000, units = "px")
print(cuadrants.plot)
dev.off()

##################################
##################################
# spatial distance and path length line plots
p_dodge <- position_dodge(0.2)

distances.data <- full.data %>% filter(path.length != 0) %>% group_by(simulation.ID,parameter.values,spatial.distance,path.length) %>% summarise(avg.net.effect = mean(abs(net.effect)))

# distances.plot <- ggplot(distances.data,aes(x = spatial.distance,y = log(avg.net.effect), group = simulation.ID))
distances.plot <- ggplot(distances.data,aes(x = spatial.distance,y = avg.net.effect, group = interaction(simulation.ID,parameter.values)))

distances.plot <- distances.plot + stat_summary(aes(color = simulation.ID),
                                                fun.y = mean,
                                                fun.ymin = function(x) ifelse(mean(x) - sd(x) > 0,mean(x) - sd(x),0), 
                                                fun.ymax = function(x) mean(x) + sd(x), 
                                                geom = "errorbar",
                                                position = p_dodge, width = 0.2) 

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

distances.plot <- distances.plot + facet_wrap(simulation.ID~parameter.values,drop = TRUE,ncol = 3, scales = "free")
distances.plot <- distances.plot + scale_color_manual(name = "", values = cbPalette[c(2,4,6)])
# distances.plot <- distances.plot + scale_color_grey(name = "",start = 0.8,end = 0.2)
# distances.plot <- distances.plot + theme_Publication() 
distances.plot <- distances.plot + xlab("spatial distance") + ylab("net effect") + guides(color = FALSE)
# distances.plot <- distances.plot + ylim(0,0.26)
# distances.plot <- distances.plot + scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05))


tiff(paste("./results/images/Fig_S31_v",file.version,".tiff",sep=""), res=600, compression = "lzw", width = 6000, height = 5000, units = "px")
print(distances.plot)
dev.off()

###############
###############

# averaged by path length

# path.length.plot <- ggplot(distances.data,aes(x = path.length,y = log(avg.net.effect), group = simulation.ID))
path.length.plot <- ggplot(distances.data,aes(x = path.length,y = avg.net.effect, group = interaction(simulation.ID,parameter.values)))

path.length.plot <- path.length.plot + stat_summary(aes(color = simulation.ID),
                                                    fun.y = mean,
                                                    fun.ymin = function(x) ifelse(mean(x) - sd(x) > 0,mean(x) - sd(x),0), 
                                                    fun.ymax = function(x) mean(x) + sd(x), 
                                                    geom = "errorbar",
                                                    position = p_dodge, width = 0.2) 

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

path.length.plot <- path.length.plot + facet_wrap(simulation.ID~parameter.values,drop = TRUE,ncol = 3, scales = "free")
path.length.plot <- path.length.plot + scale_color_manual(name = "", values = cbPalette[c(2,4,6)])
# path.length.plot <- path.length.plot + scale_color_grey(name = "",start = 0.8,end = 0.2)
# path.length.plot <- path.length.plot + theme_Publication()
path.length.plot <- path.length.plot + xlab("average path length") + ylab("net effect") + guides(color = FALSE)
# path.length.plot <- path.length.plot + scale_x_continuous(limits = c(0.85,7.2),breaks = 1:7)
# path.length.plot <- path.length.plot + scale_y_continuous(limits = c(0,0.26),breaks = seq(0,0.25,0.05))
# path.length.plot <- path.length.plot + ylim(0,0.26)

tiff(paste("./results/images/Fig_S32_v",file.version,".tiff",sep=""), res=600, compression = "lzw", width = 6000, height = 5000, units = "px")
print(path.length.plot)
dev.off()


