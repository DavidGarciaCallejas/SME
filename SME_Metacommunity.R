
#' Generate a metacommunity matrix
#'
#' @param A interaction matrix with quantitative coefficients
#' @param S number of species
#' @param N number of communities
#' @param meta.sd.a standard deviation of the interspecific coefficients among communities
#' @param meta.sd.d standard deviation of the intraspecific coefficients among communities
#' @param type either "top","bottom_nostruct", or "bottom_struct".Top: random variability in 
#' interspecific interaction coefficients among communities, with mean A and sd meta.sd.a for each coefficient.
#' bottom_nostruct: random variability in intraspecific interaction coefficients within and among communities,
#' with mean -1 and sd meta.sd.a.
#' bottom_struct: random variability in intraspecific interaction coefficients among -but not within- communities,
#' with mean -1 and sd meta.sd.a.
#' @return NxS matrix
#' @export
#'
#' @examples
Metacommunity_SME <- function(A,S,N,meta.sd.a,meta.sd.d,type){
  
  # are coefficients symmetric? check A
  symmetric.coefs <- FALSE
  if(abs(sum(A[lower.tri(A)])) == abs(sum(A[upper.tri(A)]))){
    symmetric.coefs <- TRUE
  }
  
  if(type=="top") {
    metaA = matrix(0,nr=N*S,nc=N*S)
    for(i in 1:N) {
      if(i == 1) metaA[1:S,1:S] = A
      else {
        # maintain the sign structure on each community
        # 1 - generate random from A deviates without sign, zero-truncated
        M3 = matrix(truncnorm::rtruncnorm(n = S^2,a = 0,mean = abs(A),sd = meta.sd.a),nr=S,nc=S)
        M3[A==0] = 0
        if(symmetric.coefs){
          M3[upper.tri(M3)] <- t(M3)[upper.tri(M3)]
        }
        # multiply by the sign structure of A
        M3 <- M3*sign(A)
        metaA[((i-1)*S+1):(i*S),((i-1)*S+1):(i*S)] = M3
      }			
    }	
    # Attribute -1 coefficients to the diagonal
    diag(metaA) = -1
  }
  
  else if(type=="bottom_nostruct") {
    metaA = matrix(0,nr=N*S,nc=N*S)
    for(i in 1:N) {
      if(i == 1) metaA[1:S,1:S] = A
      else {
        metaA[((i-1)*S+1):(i*S),((i-1)*S+1):(i*S)] = A
      }		
    }	
    # Attribute coefficients to the diagonal
    diag(metaA) = rnorm(N*S,-1,meta.sd.d)
  }		
  
  else if(type=="bottom_struct") {
    metaA = matrix(0,nr=N*S,nc=N*S)
    for(i in 1:N) {
      if(i == 1) metaA[1:S,1:S] = A
      else {
        metaA[((i-1)*S+1):(i*S),((i-1)*S+1):(i*S)] = A
      }		
    }	
    # Attribute coefficients to the diagonal
    E = rnorm(N,-1,meta.sd.d)
    diag(metaA) = expand.grid(c(1:S),E)[,2]
  }		
  return(metaA)
}
  
  
  