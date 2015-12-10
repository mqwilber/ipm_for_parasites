###############################################################################
## Supplement 1 for Ecological Archives: R script for calculating eviction measures
## Williams, Miller and Ellner
## Avoiding uninentional eviction from integral projection models
## (Ecology, In Press (as a Statistical Note), May 2012)
#################################################################################

# The code provided here defines functions "evictionMeasuresFC", which calculates epsilon(x) and dLambda
# as defined in the text and the Appendix, and "evictionMeasures RA", which is defined in the Appendix, as an
# alternative approach for calculating dLambda. "evictionMeasuresFC.Iter" uses iteration to calculate
# the dominant eigenvalue and eigenvector, instead of eigen(). Iteration is preferable to eigen()
# for large matrices. See Ellner and Rees (2006, Am Nat), p. 416. For smaller IPM matrices
# (250 x 250 or less) iteration is fast and eigen is faster, but for (1000 x 1000) or larger
# iteration is much faster than eigen and it can handle larger matrices.

require(Matrix);

#################################################################################
## Function to compute (lambda,w) by iteration
domeig=function(A,tol=1e-8) {
    qmax=10*tol; lam=1; x=rep(1,nrow(A));
    A=Matrix(A);
    while(qmax>tol) {
        x1=A%*%x;
        qmax=sum(abs(x1-lam*x));
        lam=sum(x1);
        x=x1/lam;
    }
    return(list(lambda=lam,w=x/sum(x)))
}

############################################################################
evictionMeasuresFC.Iter=function(growthKernel, kernel, survivalFunction,
        minsize, maxsize, parms, n.big.matrix=250) {
# Eviction measures for "floor-ceiling" solution to eviction,
# using iteration to compute dominant eigenvalue/vectors.
#
# Computes approximation to the effect of expanding the size range from
# (minsize, maxsize) by applying a demographic floor at minsize and
# demographic ceiling at maxsize. Computes everything internally, given the
# kernels & size range, using midpoint rule.
# REQUIRED ARGUMENTS
#  growthKernel = growth kernel function g(new.size, old.size, parms)
#  kernel = complete kernel function K(new.size,old.size,parms)
#  minsize,maxsize = size limits in the model, aka [L,U] or [xmin,xmax]
#  parms = parameter vector passed to the kernels. parms must be defined
#      and used in the call even if it is not used by kernel or growthKernel.
#  survivalFunction = size-dependent survival function s(x).
# OPTIONAL ARGUMENTS
#  n.big.matrix = linear dimension of approximating matrix for the kernel
############################################################################
  growthKernel=match.fun(growthKernel);
    survivalFunction=match.fun(survivalFunction);
    kernel=match.fun(kernel);

    h = (maxsize-minsize)/n.big.matrix;
    y = minsize + h*c(1:n.big.matrix)-(h/2);

    # kernel,v,w for uncorrected model
    Kmat=matrix(0,n.big.matrix,n.big.matrix);
    for(j in 1:n.big.matrix) {
        Kmat[,j]=h*kernel(y,y[j],parms);
    }
    domeigK=domeig(Kmat);
    lambda=domeigK$lambda; w=domeigK$w;
    v=domeig(t(Kmat))$w;

    eps = 1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), minsize, maxsize)$value);
    eps.U = sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), maxsize, Inf)$value);
    eps.L = sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), -Inf, minsize)$value);

    sx=survivalFunction(y,parms);
    rho=eps*sx; rho.U=eps.U*sx; rho.L=eps.L*sx;

    # Construct the iteration matrix for the corrected model
    # add big and small classes as 2 columns at the right
    Kmat2=cbind(Kmat, kernel(y,maxsize,parms),kernel(y,minsize,parms));

    # add the bottom rows: evictees are sent to the small or large class
    Kmat2=rbind(Kmat2, c(h*rho.U,rho.U[n.big.matrix],rho.U[1]));
    Kmat2=rbind(Kmat2, c(h*rho.L,rho.L[n.big.matrix],rho.L[1]));

    lambda2=domeig(Kmat2)$lambda;

    vnew=v+(1/lambda)*(v[n.big.matrix]*rho.U+v[1]*rho.L);
    dlambdaU=vnew[n.big.matrix]*sum(rho.U*w)/sum(vnew*w);
    dlambdaL=vnew[1]*sum(rho.L*w)/sum(vnew*w);

    return(list(y=y, evict=eps,evict.U=eps.U,evict.L=eps.L,lambda=lambda,
        lambda2=lambda2,dlambda=lambda2-lambda,
        dlambdaU=dlambdaU,dlambdaL=dlambdaL))
}

##############################################################################
## Eviction measures using floor - ceiling method, using eigen to calculate
## lambda, v, and w. For small matrices this is faster than using iteration
##############################################################################
evictionMeasuresFC=function(growthKernel, kernel, survivalFunction,
                            minsize, maxsize, parms, n.big.matrix=500) {
  # Eviction measures for "floor-ceiling" solution to eviction.
  #
  # Computes approximation to the effect of expanding the size range from
  # (minsize, maxsize) by applying a demographic floor at minsize and
  # demographic ceiling at maxsize. Computes everything internally, given the
  # kernels & size range
  # REQUIRED ARGUMENTS
  #  growthKernel = growth kernel function g(new.size, old.size, parms)
  #  kernel = complete kernel function K(new.size,old.size,parms)
  #  minsize,maxsize = size limits in the model, aka [L,U] or [xmin,xmax]
  #  parms = parameter vector passed to the kernels. parms must be defined
  #      and used in the call even if it is not used by kernel or growthKernel.
  #  survivalFunction = size-dependent survival function s(x).
  # OPTIONAL ARGUMENTS
  #  n.big.matrix = linear dimension of approximating matrix for the kernel
  ############################################################################
  growthKernel=match.fun(growthKernel);
  survivalFunction=match.fun(survivalFunction);
  kernel=match.fun(kernel);

  # define mesh for midpoing rule
  h = (maxsize-minsize)/n.big.matrix;
  b = minsize+c(0:n.big.matrix)*h;
  y = 0.5*(b[1:n.big.matrix]+b[2:(n.big.matrix+1)]);

  # kernel,v,w for uncorrected model
  Kmat = h*outer(y,y,kernel,parms=parms);
  w=eigen(Kmat)$vectors[,1]; w=abs(w)/sum(abs(w));
  v=eigen(t(Kmat))$vectors[,1]; v=abs(v)/max(abs(v));

  eps = 1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), minsize, maxsize)$value);
  eps.U = sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), maxsize, Inf)$value);
  eps.L = sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), -Inf, minsize)$value);

  sx=survivalFunction(y,parms);
  rho=eps*sx; rho.U=eps.U*sx; rho.L=eps.L*sx;

  # Construct the iteration matrix for the corrected model
  # add small and big classes as 2 columns at the right
  Kmat2=cbind(Kmat, kernel(y,minsize,parms),kernel(y,maxsize,parms));

  # add the bottom rows: evictees are sent to the small or large class
  Kmat2=rbind(Kmat2, c(h*rho.L,rho.L[1],rho.L[n.big.matrix]));
  Kmat2=rbind(Kmat2, c(h*rho.U,rho.U[1],rho.U[n.big.matrix]));
  lambda=abs(eigen(Kmat)$values[1]);
  lambda2=abs(eigen(Kmat2)$values[1]);

  vnew=v+(1/lambda)*(v[n.big.matrix]*rho.U+v[1]*rho.L);
  dlambdaU=vnew[n.big.matrix]*sum(rho.U*w)/sum(vnew*w);
  dlambdaL=vnew[1]*sum(rho.L*w)/sum(vnew*w);

  return(list(y = y, evict=eps,evict.U=eps.U,evict.L=eps.L,lambda=lambda,
              lambda2=lambda2,dlambda=lambda2-lambda,
              dlambdaU=dlambdaU,dlambdaL=dlambdaL))
}

###############################################################################
## EVICTION MEASURES using the reassignment method (Solution 2 in Appendix)
###############################################################################
evictionMeasuresRA=function(growthKernel, kernel, minsize, maxsize, parms,
                    n.big.matrix=250, minfate=NULL, maxfate=NULL, survivalFunction) {
  # REQUIRED ARGUMENTS
  #  growthKernel = growth kernel function g(new.size, old.size, parms)
  #  kernel = complete kernel function K(new.size,old.size,parms)
  #  minsize,maxsize = size limits in the model, aka [L,U] or [xmin,xmax]
  #  parms = parameter vector passed to the kernels. parms must be defined
  #      and used in the call even if it is not used by kernel or growthKernel,
  #      in the same way that parms must be used in lsoda().
  # OPTIONAL ARGUMENTS
  #  n.big.matrix = linear dimension of approximating matrix for the kernel,
  #     i.e., number of size intervals used for integration by midpoint rule
  #  (minfate,maxfate) = range of plausible future sizes for evictees;
  #    these default to (minsize, maxsize)
  #  survivalFunction = size-dependent survival function s(x). If not provided,
  #    $d\lambda$ is not computed
  ############################################################################
  growthKernel=match.fun(growthKernel);
  survivalFunction=match.fun(survivalFunction);
  kernel=match.fun(kernel);
  h = (maxsize-minsize)/n.big.matrix;
  b = minsize+c(0:n.big.matrix)*h;
  y = 0.5*(b[1:n.big.matrix]+b[2:(n.big.matrix+1)]);
  Kmat = h*outer(y,y,kernel,parms=parms);
  eK=eigen(Kmat);
  w=eK$vectors[,1]; w=abs(w)/sum(abs(w));
  lambda=abs(eK$values[1]);
  v=eigen(t(Kmat))$vectors[,1]; v=abs(v)/max(abs(v));
  eps=1-sapply(y,function(z) integrate(function(u) growthKernel(u,z,parms), minsize, maxsize)$value);
  if(is.null(minfate)) minfate=minsize;
  if(is.null(maxfate)) maxfate=maxsize;
  e=(y<minfate)|(y>maxfate)
  vmax=max(v[!e]);
  sx=survivalFunction(y,parms);
  vnew=v+sx*eps*vmax/lambda;
  vmax=max(vnew[!e]);
  dL1=vmax*sum(sx*eps*w)/sum(vnew*w);
  return(list(evict=eps,dLambda=dL1,lambda=lambda))
}