#' Generate an interaction matrix with several constraints
#' 
#' It accepts a one-dimensional niche axis and a mode that specifies what niche positions are more prone to interact
#' Interactions can be more likely to happen between close-by species (mode = "same"), 
#' species in adjacent positions (mode = "adjacent),
#' or distant species (mode = "other")
#'
#' @param S number of species
#' @param connectance connectance of the matrix
#' @param niches one-dimensional niche axis, a numeric vector of length S
#' @param mode which niche position is more likely to interact: same, adjacent, other
#' @param forbidden.positions either 0 or a matrix with 1s in the forbidden interactions
#' @param connectance.type directed or undirected
#'
#' @return qualitative interaction matrix
#' @export
#'
#' @examples
InteractionTopology_SME <- function(S,connectance,niches,mode,forbidden.positions = 0,connectance.type = "undirected"){
  
  # keep track of the species index
  if(is.null(names(niches))){
    names(niches) <- 1:length(niches)
  }
  
  # create matrix
  interaction.matrix <- matrix(0,nrow = S,ncol = S)
  
  # how many links?
  link.number <- ifelse(connectance.type == "undirected",round(connectance*((S*(S-1))/2)),round(connectance*((S^2)-S)))
  
  if(link.number>0){
    # generate the trophic position boundaries for each species
    # and store them in a list
    lower.range <- list()
    higher.range <- list()
    lower.same <- list()
    lower.adjacent <- list()
    lower.other <- list()
    higher.same <- list()
    higher.adjacent <- list()
    higher.other <- list()
    
    for(i.sp in 1:S){
      # basically, divide the niche axis in lower and higher parts from the species' niche
      # each part is divided in 3 equal segments, that correspond 
      # to the "same","adjacent", and "other" trophic levels
      
      # then, the pdf is a normal (or a mixture thereof) 
      # centered in the segment with highest probability
      # sd is half the length of the segment
      
      lower.range[[i.sp]] <- c(0,niches[i.sp])
      higher.range[[i.sp]] <- c(niches[i.sp],1)
      
      lower.same[[i.sp]] <- c(2*((lower.range[[i.sp]][2]-lower.range[[i.sp]][1])/3),niches[i.sp])
      lower.adjacent[[i.sp]] <- c(1*((lower.range[[i.sp]][2]-lower.range[[i.sp]][1])/3),2*((lower.range[[i.sp]][2]-lower.range[[i.sp]][1])/3))
      lower.other[[i.sp]] <- c(0,1*((lower.range[[i.sp]][2]-lower.range[[i.sp]][1])/3))
      
      higher.same[[i.sp]] <- c(niches[i.sp],niches[i.sp] + 1*((higher.range[[i.sp]][2]-higher.range[[i.sp]][1])/3))
      higher.adjacent[[i.sp]] <- c(niches[i.sp] + 1*((higher.range[[i.sp]][2]-higher.range[[i.sp]][1])/3), niches[i.sp] + 2*((higher.range[[i.sp]][2]-higher.range[[i.sp]][1])/3))
      higher.other[[i.sp]] <- c(niches[i.sp] + 2*((higher.range[[i.sp]][2]-higher.range[[i.sp]][1])/3),1)
    }
    
    # generate each link
    for(i.link in 1:link.number){
      
      # 1 - pick a species randomly
      affected.sp <- sample.int(n = S,size = 1)
      
      # while loop for avoiding replicated links
      existing <- TRUE
      safety.count <- 0
      while(existing & safety.count < 1e4){
        
        # 2 - generate interactions according to the mode
        
        # 2.1 - if the mode is interactions in the same trophic position,
        # the pdf is a normal centered on the affected sp's niche with a sd of 
        # the range of the "same" range
        # 2.2 - otherwise, the pdf is a mixture normal distribution, with
        # the first component in the lower part (either "adjacent" or "other")
        # and the second component in the higher part
        # both can be chosen with equal probability
        
        if(mode == "same"){
          link.position <- truncnorm::rtruncnorm(n = 1,a = 0,b = 1,mean = niches[affected.sp],sd = (higher.same[[affected.sp]][2]-lower.same[[affected.sp]][1])/2)
        }else if(mode == "adjacent"){
          means <- c(mean(lower.adjacent[[affected.sp]]),mean(higher.adjacent[[affected.sp]]))
          sds <- c((lower.adjacent[[affected.sp]][2] - lower.adjacent[[affected.sp]][1])/2,
                   (higher.adjacent[[affected.sp]][2] - higher.adjacent[[affected.sp]][1])/2)
          
          # sample from lower or higher component?
          component <- round(runif(1,1,2))
          link.position <- rnorm(n=1,mean=means[component],sd=sds[component])
        }else{
          means <- c(mean(lower.other[[affected.sp]]),mean(higher.other[[affected.sp]]))
          sds <- c((lower.other[[affected.sp]][2] - lower.other[[affected.sp]][1])/2,
                   (higher.other[[affected.sp]][2] - higher.other[[affected.sp]][1])/2)
          
          # sample from lower or higher component?
          component <- round(runif(1,1,2))
          link.position <- rnorm(n=1,mean=means[component],sd=sds[component])
        }
        
        # 3 - which species is closer to the selected link position?
        coupled.sp <- as.integer(names(which.min(abs(niches[-affected.sp] - link.position))))
        
        # 4 - fill the matrix positions, taking into account the forbidden.positions argument
        # note that I fill both positions, regardless of link type
        if(interaction.matrix[affected.sp,coupled.sp] == 0 & 
           interaction.matrix[coupled.sp,affected.sp] == 0){
          # if there are forbidden.positions, check them
          if(sum(forbidden.positions == 1)>0){
            if(forbidden.positions[affected.sp,coupled.sp] == 0 &
               forbidden.positions[coupled.sp,affected.sp] == 0){
              existing <- FALSE
              interaction.matrix[affected.sp,coupled.sp] <- 1
              interaction.matrix[coupled.sp,affected.sp] <- 1
            }#if
            # else, no forbidden.positions
          }else{
            existing <- FALSE
            interaction.matrix[affected.sp,coupled.sp] <- 1
            interaction.matrix[coupled.sp,affected.sp] <- 1
          }
        }# if link did not exist, fill it
        
        # safety.count just to avoid infinite loops
        safety.count <- safety.count + 1
      }# while loop
      
      
    }# for each link
    
  }# if link.number>0
  return(interaction.matrix)
}


