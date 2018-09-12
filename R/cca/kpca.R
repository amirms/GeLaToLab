
center.kernel.matrix = function(km) {
  m <- nrow(km)
  kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2
  kc
}


## Matrix Interface
kpca = function(km,  features = 0, th = 1e-4, na.action = na.omit, ...)
          {
            km <- as.matrix(km)
            m <- nrow(km)
           
            ## center kernel matrix
            kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2
            
            ## compute eigenvectors
            res <- eigen(kc/m,symmetric=TRUE)
            
            if(features == 0)
              features <- sum(res$values > th)
            else 
              if(res$values[features] < th)
                warning(paste("eigenvalues of the kernel matrix are below threshold!"))
            
            pcv <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
            eig <- res$values[1:features]
            names(eig) <- paste("Comp.", 1:features, sep = "")
            rotated <- kc %*% pcv
            kcall <- match.call()
            kernelf <- km

            return(list(pcv=pcv, eig=eig, rotated=rotated, kcall=kcall, kernelf=kernelf))
          }
