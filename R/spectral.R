#Input: an adjacency matrix
#Output: a normalized laplacian matrix
#Corresponds to NCUT(A,k)
normalized.symmetric.laplacian<- function(A) {
  D = diag(rowSums(A))
  sqrt.inv.D <- solve(sqrt(D))
  L = diag(dim(A)[1])-(sqrt.inv.D%*%A%*%sqrt.inv.D)
  
  dimnames(L) <- dimnames(A)  
  L
}


normalized.random.walk.laplacian<- function(A) {
  D = diag(rowSums(A))
  inv.D <- solve(D)
  L = diag(dim(A)[1])-(inv.D%*%A)
  
  dimnames(L) <- dimnames(A)  
  L  
}


#Input: an adjacency matrix
#Output: an unnormalized laplacian matrix
#Corresponds to RATIOCUT(A,k)
unnormalized.laplacian <- function(A) {
  D = diag(rowSums(A))
  L = D-A
  
  dimnames(L) <- dimnames(A)  
  L
}

#Input: an asymmetric adjacency matrix
MQ.laplacian <- function(A) {
  
  for (i in 1:dim(A)[1])
    A[i,i] = 0
  
  rows = rowSums(A)
  cols = colSums(A)

  
  #Condition: there are only zeroes on diagonal
  D <- diag(cols + rows)
  
  L = D - A
  
  return(L)
  
  inv.D <- solve(D)
  return((inv.D%*%L))  
    
}


#Input: a laplacian matrix
#Output: k, where k-1th eigen value is small
compute.eigen.gap <- function(L) {
  y<-eigen(L, only.values= TRUE)
  y$values <- rev(y$values)
  
  
  largest_gap=0
  k = 1
  n = length(y$values)
  
  for (i in 2:(n-1) ) {
    gap = y$values[i+1] - y$values[i]
    
    if (gap > largest_gap)
    {
      largest_gap <- gap
      
      k <- i
    }
    
  }
   
  return(k)
  
}


forward.derivative <- function(vec, order, h=1) {
  n = length(vec)
  
  if (order*h >= n)
    stop("undefined!")
  
  if (order >= n)
    stop("the order of differential does not match the no of data")
  
  diff = rep(0, n)
  
  for (i in 1:(n-order* h)) {
    
    for (j in 0:order) {
      
      sign = (-1) ^ j
      
      comb = choose(order, j)
      
      y = vec[i + (order-j)*h]
      
      diff[i] = diff[i] + sign * comb * y
      
      
    }
   diff[i] = diff[i] / (h ^2)
    
    
  }
 
  return(diff)
  
}

#order has to be even

central.derivative <- function(vec, order, h=1) {
  n = length(vec)
  
  if (order*h >= n)
    stop("undefined!")
  
  if ((order * h) %% 2 != 0)
    stop("the order of differential has to be even")
  
  diff = rep(0, n)
  
  for (i in ((h*order/2)+1):(n - (order*h/2)))
    
    for (j in 0:order) {
      
      print(i)
      
      sign = (-1) ^ j
      
      comb = choose(order, j)
      
      y = vec[i + ((order/2)-j)*h]
      
      diff[i] = diff[i] + sign * comb * y
      
      
    }
  
  return(diff)
  
}



#Input: a single dim - nrow=ncol
# prob of generating a 1.
#postcond: a square matrix
random.binary.nonsymmetric.matrix <- function(ndim, prob=0.5) {
  
  m <- apply(matrix(nrow=ndim, ncol=ndim), c(1,2), function(x) sample(c(0,1),1,
                                                                      prob=c((1-prob), prob)))
  
  return(m)
}

#Input: a single dim - nrow=ncol
#postcond: a square matrix
random.binary.symmetric.matrix <- function(ndim, prob=0.5) {
  
  m <- apply(matrix(nrow=ndim, ncol=ndim), c(1,2), function(x) sample(c(0,1),1,
                                                                      prob=c((1-prob), prob)))
  ind <- lower.tri(m)
  m[ind] <- t(m)[ind]
  
  return(m)
}

#Input: a single dim - nrow=ncol
#postcond: a square matrix
random.symmetric.matrix <- function(ndim) {
  
  m <- apply(matrix(nrow=ndim, ncol=ndim), c(1,2), function(x) sample(0:(ndim-1),1))
  ind <- lower.tri(m)
  m[ind] <- t(m)[ind]
  
  return(m)
}


test.eigen.gap1 <- function(ndim) {
  
  A <- random.binary.symmetric.matrix(ndim)
  L1 <- normalized.laplacian(A)
  L2 <- unnormalized.laplacian(A)
  #compute.eigen.gap(L1)
  compute.eigen.gap(L2)
  
}

norm_vec <- function(x) sqrt(sum(x^2))

# 5 iterations to find the best kmeans

spectral.clustering <- function(L, k, iter = 10) {
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  V$vectors <- V$vectors[, rev(seq_len(ncol(V$vectors)))]
  
  n <- length(V$values)

  
  
  #plot(V$vectors)
  
  #FIXME normalization has to be there for symmetrix laplacian
  
  for(i in 1:n) {
    row = V$vectors[i,]
    norm.row <- norm_vec(row)
    
    for(j in 1:k)
      V$vectors[i,j] <- V$vectors[i,j]/norm.row
    
  }
  
  
  V$vectors <- V$vectors[,-(k+1:n)]

  best_cluster = list(partition=NULL, BSS_TSS=0)

  for (i in 1:iter) {
  
    cl = kmeans(V$vectors,k, iter.max = 300, nstart = 200)

    current_BSS_TSS = cl$betweenss / cl$tot.withinss

    if (current_BSS_TSS > best_cluster$BSS_TSS) {
      best_cluster$BSS_TSS = current_BSS_TSS
      
      best_cluster$partition = cl$cluster
      
      print("changed")
          
    }   
    
  }

  return(best_cluster$partition)

}

find.eigen.gap <- function(m) {
  
  
  #can change
  #L <- MQ.laplacian(m)
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  #get rid of the smallest eigenvalue
  V$values <- V$values[-1]
  plot(V$values)
  
  
  n <- dim(m)[1]
  
  fd=matrix(0, nrow= n-2, ncol = n-1)
  
  colnames(fd) = 2:n
  
  for (i in 1:(n-2))
    fd[i,] = forward.derivative(V$values, i)

  
  col = 1
  
  for (row in (n-2):2) {
    
    if (((row-1) %% 2) == 0) { 
      #even derivative, smaller is better
      
      
      if (fd[row, col] > 0) 
        col <- col+1 #go to right
      #TODO what about if it was zero
      else
       # if (abs(fd[row-1, col+1]) > abs(fd[row-1, col]) )
        if (fd[row-1, col+1] > fd[row-1, col])
          col <- col+1
      #else col remains the same
        
    }else{   
      #odd derivative, greater is better 
      if (fd[row, col] < 0) 
        col <- col+1 #go to right
      # else  go to left => col remains the same
      else
        #if (abs(fd[row-1, col+1]) < abs(fd[row-1, col]) )
        if (fd[row-1, col+1] < fd[row-1, col])
          col <- col+1
      #else col remains the same
    
    }
    
    
  }
  
  
  return(col+1)
  
}

best.fit.loess <- function(m)
{
  #can change
  #L <- MQ.laplacian(m)
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  #get rid of the smallest eigenvalue
  #get rid of the zero eigenvalues
  y <- V$values[which(V$values > 0.00001)]
  
  y = y[1:ceiling(y/2)]
  x = (length(V$values) - length(y) + 1): length(V$values)

  print(length(y)==length(x))
  
  plot(x,y,type="l",ylim=c(0,2))
  lo <- loess(y~x,span=0.05)
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  out = predict(lo,xl)
  lines(xl, out, col='red', lwd=2)
  
  infl <- c(FALSE, diff(diff(out)>0)!=0)
  
  print(infl)
  
  points(xl[infl ], out[infl ], col="blue")
  
  return(out)
  
}

best.fit.curve <- function(m) {
  #can change
  #L <- MQ.laplacian(m)
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  #get rid of the zero eigenvalues
  y <- V$values[which(V$values > 0.00001)]
  
  y = y[1:ceiling(y/2)]
  x = (length(V$values) - length(y) + 1): length(V$values)
  
  
  #y = c(1, 1.5, 1.75, 2, 2, 3, 5, 8, 12,17,23, 30, 38, 47)
  #x = 1:length(y)
  
  #third degree
  fit <- lm(y~poly(x,3,raw=TRUE))
  
  print(fit)
  
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  
  plot(x,y,type="l",ylim=c(min(y),max(y)))
  #plot(x,y,type="l",ylim=c(0,2))
  
  lines(xl, predict(fit, data.frame(x=xl)), col='red', lwd=2)
  
  #fit(x) = ax^3 + bx^2 + cx + d
  
  a = fit$coefficients[4]
  b = fit$coefficients[3]
  c = fit$coefficients[2]
  d = fit$coefficients[1]
  
  inflection = ((-1) * b)/(3*a)
  
  discrete.inflection = floor(inflection)
  
  print("inflection point is:")
  print(inflection)
  
  residuals = fit$residuals
    
  pos.neg = 1
  neg.pos = 1
  
  for (i in (discrete.inflection-1):2)
    if ((residuals[i] > 0) && (residuals[i-1]<0)) {
      pos.neg = x[i-1]
      break     
    }
  
  print(pos.neg)
  

  for (i in (discrete.inflection-1):2)
    if ((residuals[i] < 0) && (residuals[i-1]>0)) {
      neg.pos = x[i-1]
      break      
    }
  
  print(neg.pos)
  
  pos.neg.mq = compute.MQ(m, pos.neg)  
  print(pos.neg.mq$MQ)
  
  neg.pos.mq = compute.MQ(m, neg.pos)
  print(neg.pos.mq$MQ)  
  
  return(residuals)
  
}



best.fit.quadratic.curve <- function(m) {
  #can change
  #L <- MQ.laplacian(m)
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  #get rid of the zero eigenvalues
  zeroes = which(V$values < 0.00001)
  
  y = V$values
  
  x = (length(zeroes)+1):ceiling(length(y)/2)
  
  print(x)
  
  y = y[(length(zeroes)+1):ceiling(length(y)/2)]
  
  #y = c(1, 1.5, 1.75, 2, 2, 3, 5, 8, 12,17,23, 30, 38, 47)
  #x = 1:length(y)
  
  #third degree
  fit <- lm(y~poly(x,2,raw=TRUE))
  
  print(fit)
  
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  
  plot(x,y,type="l",ylim=c(min(y),max(y)))
  #plot(x,y,type="l",ylim=c(0,2))
  
  lines(xl, predict(fit, data.frame(x=xl)), col='red', lwd=2)
  
  #fit(x) = ax^3 + bx^2 + cx + d
  
  #a = fit$coefficients[4]
  #b = fit$coefficients[3]
  #c = fit$coefficients[2]
  #d = fit$coefficients[1]
  
  #inflection = ((-1) * b)/(3*a)
  
  #discrete.inflection = floor(inflection)
  
  
  return(fit)
  
}


fit.spectral.embedding <- function(m) {
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  V<-eigen(L)
  #V$values, V$vectors
  k=2
  V$values <- rev(V$values)
  
  V$vectors <- V$vectors[, rev(seq_len(ncol(V$vectors)))]
  
  n <- length(V$values)

  for(i in 1:n) {
    row = V$vectors[i,]
    norm.row <- norm_vec(row)
    
    for(j in 1:k)
      V$vectors[i,j] <- V$vectors[i,j]/norm.row
    
  }
  
  return(V$vectors[,1:2])
  
}


derivatives.vs.MQ <- function(m, order=3) {
  
  require(scatterplot3d)
  
  #can change
  #L <- MQ.laplacian(m)
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  #get rid of the smallest eigenvalue
  V$values <- V$values[-1]
  
  dev = V$values
  
  n <- dim(m)[1]
  
  fd=matrix(0, nrow= n-1, ncol = n-1)
  cd=matrix(0, nrow = (n+1)/2, ncol = n-1)
  #cd=matrix(0, nrow= n-1, ncol = n-1)
  
  colnames(fd) = 2:n
  fd[1,] = dev
  
  colnames(cd) = 2:n
  cd[1,] = dev
  
  #cd[2,] = forward.derivative(dev, 1)
  
  print(n)
  
  for (i in 1:((n-1)/2)) {

    #fd[i+1,] = forward.derivative(dev, i)
  
    #if ((i %% 2) == 0)  {
      
    print(i)
    
    v= central.derivative(dev, 2, i)
    
    print(v)
    
      cd[i+1,] = v
      
    #}
      
  } 
  
  return(cd)
  
  
  #curvature = fd[2,] / fd[3,]
  
  #indices = which(curvature < 0)
  
  
  return(fd)
    
  scores = rep(0, length(indices))
  
  for (i in 1:length(indices)) {
    
    #print(indices[i])
    
    partition <- spectral.clustering(L, indices[i]+1)
    
    names(partition) <- colnames(m)
    
    partition <- normalizeVector(partition)
    
    scores[i] = evaluateMQ(m, partition) 
  }
  
  print(fd[5, indices])
  
  scatterplot3d(fd[12, indices], indices, scores,        # x y and z axis
                       color="blue", pch=10,        # filled blue circles
                       type="h",                    # vertical lines to the x-y plane
                       xlab="no of clusters",
                       ylab="11th Derivative",
                       zlab="MQ Scores")
  
  #plot(fd[4, indices], scores)

}


#Input: a list of k's: K

compute.MQ <- function(m, K) {
  #can change
  #L <- MQ.laplacian(m)
  #L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  #REMOVE ME LATER
  V<-eigen(L)
  #V$values, V$vectors
  
  V$values <- rev(V$values)
  
  plot(V$values)
  
  n <- dim(m)[1]
  
  best_cluster = list(partition=NULL, k =1, MQ=0)
  
  for (i in 1:length(K)) {
    
    partition <- spectral.clustering(L, K[i])
    
    #print("Total within-cluster sum of squares")
    #print(P$tot.withinss)
    
    #print("The between-cluster sum of squares")
    #print(P$betweenss)
    
    #print("BSS/TSS")
    
    #print(P$betweenss / P$tot.withinss)
    
    names(partition) <- colnames(m)
    
    partition <- normalizeVector(partition)
    
    current_k = max(partition)
    
    print(current_k)
    
    current_MQ = evaluateMQ(m, partition) 
    
    print(current_MQ)
    
    
    if (current_MQ > best_cluster$MQ) {
      best_cluster$MQ = current_MQ
      
      best_cluster$partition = partition
      
      best_cluster$k = current_k
      
      
    }   
    
    
  }
 
  
  return(best_cluster)
  
  
}


compute.random.MQ <- function(m, size) {
  sm <- make.symmetric(m)
  
  n <- dim(m)[1]
  
  s <- sample(1:n, size)
  
  #can change
  L <- normalized.symmetric.laplacian(sm)
  
  for (i in 1:size) {
    partition <- spectral.clustering(L, s[i])
    
    names(partition) <- colnames(m)
    
    partition <- normalizeVector(partition)
    
    print(max(partition))
    
    print(evaluateMQ(m, partition))
     
    
  }
  
}


bisect <- function(partition, indices) {
  k = max(partition)
  
  for (i in 1:k) {
    
    i.indices = which(partition==i)
    
    partition[which(i.indices %in% indices)]= k+i
    
  }
  
  partition <- normalizeVector(partition)
  
  return(partition)
  
}

recursive.bisection <- function(m) {
  
  n <- dim(m)[1]
  
  best = list(partition = rep(1,n), MQ=0)
  
  best$MQ = evaluateMQ(m, best$partition)
  
  #can change
  L <- MQ.laplacian(m)
  #L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(m)
  
  V <- eigen(L)
  
  V$values <- rev(V$values)
  
  V$vectors <- V$vectors[, rev(seq_len(ncol(V$vectors)))]
  
  #sort ascendingly
  
  for (i in 2:n) {
    
    u = V$vectors[,i]
    
    t = median(u)

    indices = which(u > t)
    
    p <- bisect(best$partition, indices)
    
    current_MQ <- evaluateMQ(m, p)
    
    print(current_MQ)
    
    if(best$MQ > current_MQ)
      break
    
    best$partition = p
    best$MQ = current_MQ
    
  }
  
  
  return(best)
  
}

make.symmetric <- function(m) {
  pmean <- function(x,y) (x+y)/2
  
  m[] <- pmean(m, matrix(m, nrow(m), byrow=TRUE))
  
  return(m)
}


best.fit <- function(m) {
  
  L <- normalized.symmetric.laplacian(make.symmetric(m))
  #L <- unnormalized.laplacian(make.symmetric(m))
  #L <- normalized.random.walk.laplacian(make.symmetric(m))
  
  n <- dim(L)[1]
  
  V <- eigen(L)
  
  
  
  V$values <- rev(V$values)
  
  V$vectors <- V$vectors[, rev(seq_len(ncol(V$vectors)))]
  
  y = V$values
  x = 1:length(V$values)
  
  plot(x,y, pch=5, xlim=c(0,n))
  
  fit <- lm(y~poly(x,50,raw=TRUE))
  
#  xx <- seq(0,n, length=50)
  
  lines(x, predict(fit, data.frame(x)), col="red")

  return(fit)
}

compute.prob <- function(A) {
  #outgoing edges
  D = rowSums(A)
  
  for(i in 1:dim(A)[1])
    for(j in 1:dim(A)[2])
      if (A[i,j] > 0)
        A[i,j] = A[i,j] / D[i]
  
  return(A)
}

compute.fix.point <- function(x0, A, alpha= 1) {
  
  P = compute.prob(A)
  
  y = eigen(P)
  
  return(P)
  
  print(y)
  
  return(x0 %*% solve(iden - alpha * A))
  
}

