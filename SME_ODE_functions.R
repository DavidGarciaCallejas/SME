# Function to evaluate the steady-state
model_stode = function(t,y,parms=NULL,A,D) {
	dy = y*(1+A%*%y) + D%*%y
	return(list(dy,1))
}

# Function to evaluate the Jacobian matrix	
model_J = function(t,y,parms=NULL,A,D) {
	dy = y*(1+A%*%y) + D%*%y
	return(as.list(dy))
}

# Function for progressive finding of the solution for local communities
eq_local = function(metaA,S,N,nstep) {	
	
	SSmat = matrix(nr=N,nc=S)

	for(x in 1:N) {
		if(x == 1) A = metaA[1:S,1:S] 
		else A = metaA[((x-1)*S+1):(x*S),((x-1)*S+1):(x*S)]
						
		# Steady state with only intraspecific competition
		A0 = diag(diag(A))
		SSmat[x,] = -diag(A0)^-1	

		# Progressive solution
		for(i in 0:nstep) {
			A1 = A0 + A*i/nstep 
			diag(A1) = diag(A)
			SSmat[x,] = stode(y = SSmat[x,],time=0,func=model_stode,parms=NULL,A=A1,D=matrix(0,S,S), positive = TRUE)[[1]] 
		}		

	}
	
	SSvec = NULL
	for(x in 1:N) SSvec = c(SSvec,SSmat[x,])
	
	return(SSvec)
}










