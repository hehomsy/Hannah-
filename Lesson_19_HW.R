options(warn=-1)
library(pracma)
#### Outer  - create contour plots ####################################
##########################################################################

Outer = function(f,x){
  n1 = length(x[[1]])
  n2 = length(x[[2]])
  res = matrix(0,nrow=n1,ncol=n2)
  rownames(res) = x[[1]]
  colnames(res) = x[[2]]
  for (i in 1:n1){
    for (j in 1:n2){
      res[i,j]=f(c(x[[1]][[i]],x[[2]][[j]]))
    }
  }
  return(res)
}


##########################################################################
#### meshgrid  - used by vector field ####################################
##########################################################################
meshgrid = function(x,y){
  n = length(x)
  m = length(y)
  mat1 = matrix(0,nrow=m,ncol=n)
  mat2 = matrix(0,nrow=m,ncol=n)
  for (i in 1:m){
    mat1[i,]=x
  }
  for (i in 1:n){
    mat2[,i]=y
  }
  ans = list(X=mat1,Y=mat2)
  return (ans)
}

##########################################################################
#### Quiver  - used by vector field ######################################
##########################################################################

quiver <- function(x, y, u, v,
                   scale = 0.05, angle = 10, length = 0.1, ...) {
  stopifnot(is.numeric(x), is.numeric(y), is.numeric(u), is.numeric(v))
  
  arrows(x, y, x+scale*u, y+scale*v, angle=10, length=length, ...)
}

###########################################################################
##### Vector Field ####################
###########################################################################

VectorField = function(fun, xlim, ylim, n = 16,
                       scale = 0.05, col = "darkblue",xlab = "xlim", ylab="ylim",
                       main="",...) {
  stopifnot(is.numeric(xlim), length(xlim) == 2,
            is.numeric(ylim), length(ylim) == 2)
  
  xpts = seq(xlim[1],xlim[2],length.out=n)
  ypts = seq(ylim[1],ylim[2],length.out=n)
  
  M = meshgrid(xpts, ypts)
  
  x = M$X
  y = M$Y
  px=M$X
  py=M$Y
  for (i  in 1:n){
    for (j in 1:n){
      ans = fun(c(xpts[j],ypts[i]))
      px[i,j]=ans[1]
      py[i,j]=ans[2]
    }
  }
  
  
  plot(xlim, ylim, type="n",xlab=xlab,ylab=ylab,main=main); grid()
  quiver(x, y, px, py, scale = scale, col = col, ...)
  #return(list(px=px,py=py))
}

###########################################################################
##### Jacobian of a 2D vector function ####################################
###########################################################################
Jacobian2 = function(f,x0,h=1E-4){
  jax = matrix(0,nrow=2,ncol=2)
  xph = c(x0[1]+h,x0[2]);xmh=c(x0[1]-h,x0[2])
  yph = c(x0[1],x0[2]+h);ymh=c(x0[1],x0[2]-h)
  jax[,1]=(f(xph)-f(xmh))/(2*h)
  jax[,2]=(f(yph)-f(ymh))/(2*h)
  return(jax)
}

###########################################################################
##### Zeros - root finding for a vector function ##########################
###########################################################################
Zeros = function(f,x0,h=1E-4,tol=1E-4){
  i = 1
  p=c(1,1)
  while (Norm(p)>tol & i<100){
    p = solve(Jacobian2(f,x0),-f(x0)) # linear algebra step
    x0 = x0+p
    #print(i)
    #print(x0)
    i=i+1
  }
  return(x0)
}

###########################################################################
##### path - output is the points along a particular dynamic system########
###########################################################################
path = function(f,x0,deltat=0.01,N=1000,tol=1E-4){
  len = length(x0)
  points=matrix(0,ncol=len)
  points[1,] = x0
  n = 0
  p = c(1,1)
  while(Norm(p)>tol & n<N){
    n=n+1
    p = f(x0)*deltat
    x0=x0+p
    points = rbind(points,x0)
  }
  
  rownames(points)=0:n
  return(points)
}

##Problem 1

#a lamba = 3 w =.25, 
f = function(x){c((-w*lambda*0.05*x[2]-lambda*(0.005)*x[1]*x[2])*((x[1]>0)*(x[2]>0)),
                  (-w*0.05*x[1]-0.005*x[1]*x[2])*((x[1]>0)*(x[2]>0)))}
lambda = 3
w = .25

f = function(lambda,w){
  f1 = function(x){c((-w*lambda*0.05*x[2]-lambda*(0.005)*x[1]*x[2])*((x[1]>0)*(x[2]>0)),
                    (-w*0.05*x[1]-0.005*x[1]*x[2])*((x[1]>0)*(x[2]>0)))}
  x=c(5,2)
  pts = path(f1,x,deltat=1)
  N = length(pts[,1])
  if (pts[N,1]>pts[N,2]){
    winner = "R"
    hours = min(which(pts[,2]<0))-1
    forces = as.numeric(round(pts[N,1],1))#round(pts[N,1],1)
  }
  else{
    winner="B"
    hours = min(which(pts[,1]<0))-1
    forces = as.numeric(round(pts[N,2],1))
  }
  
  return (list(winner=winner,hours=hours,forces=forces))
  
}
ans = f(3,.25)
print(ans)  #Red forces win and that takes 39 hours, they will have 2.2 forces left

#b repeat for w = .1, .2, .5, .75, .9
lambda = 3
w = seq(.1,.9,.05)
ans.hours=0
ans.win = "R"
ans.forces = 0
for (i in 1:length(w)){
  res = f(lambda,w[i])
  ans.hours[i]=res$hours
  ans.win[i]=res$winner
  ans.forces[i]=res$forces 
}
result=data.frame(w=w,hours=ans.hours,winner=ans.win,forces=ans.forces)
print(result)

#c R side benefits from fighting in adverse weather conditions, so if I was the blue commander, I would expect red to attack on a rainy day

#d
lambda = seq(1.5,5,.5)
w = .25
ans.hours=0
ans.win = "R"
ans.forces = 0
for (i in 1:length(lambda)){
  res = f(lambda[i],w)
  ans.hours[i]=res$hours
  ans.win[i]=res$winner
  ans.forces[i]=res$forces 
}
result=data.frame(lambda=lambda,hours=ans.hours,winner=ans.win,forces=ans.forces)
print(result)

lambda = seq(1.5,5,.5)
w = w = seq(.1,.9,.05)
ans.hours=0
ans.win = "R"
ans.forces = 0
for (i in 1:length(lambda)){
  for(j in 1:length(w)){
    res = f(lambda[i],w[j])
    ans.hours[i][j]=res$hours
    ans.win[i][j]=res$winner
    ans.forces[i][j]=res$forces 
  }
}
result=data.frame(lambda=lambda,w=w,hours=ans.hours,winner=ans.win,forces=ans.forces)
print(result)

