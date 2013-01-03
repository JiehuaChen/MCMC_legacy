incor <- function(d.site, phi.draw, sph.cor20.site){
	#stack diagonally
	d.site.scaled <- lapply(lapply(d.site, "-"), "/", phi.draw) 
	cor <- lapply(d.site.scaled, exp)
	cor.taper <- mapply("*", cor, sph.cor20.site)
	invcor <- lapply(cor.taper, solve)
	return(invcor)
}

cov.full <- function(d.site, phi.draw, sph.cor20.site,  sigma2.alpha.draw, sigma2.draw,sizes.noise, sizes.site, size.points){
	d.site.scaled <- lapply(lapply(d.site, "-"), "/", phi.draw) 
	cor <- lapply(d.site.scaled, exp)
#stack diagonally
	sizes.site <- lapply(d.site, dim)
	cov.stack <- matrix(0, size.points, size.points)
	cov.stack[(1:sizes.noise), (1:sizes.noise)] <- (sigma2.draw+sigma2.alpha.draw)*diag(rep(1, sizes.noise))
	k <- 1+sizes.noise
	for(i in 1:length(d.site)){
		cov.stack[(k:(k+sizes.site[[i]][1]-1)), (k:(k+sizes.site[[i]][1]-1))] <- c(sigma2.alpha.draw)*cor[[i]]+c(sigma2.draw)*diag(rep(1, sizes.site[[i]][1]))
		k <- k+(sizes.site[[i]])[1]
	}
	return(cov.stack)
}

quardterm <- function(d.site, alpha.draw,sizes.noise, inv.cor.draw){
	sizes.site <- lapply(d.site, dim)
	termvalue <- 0
	k <- sizes.noise+1
	for(i in 1:length(d.site)){
		termvalue <- termvalue+ t(alpha.draw[k:(k+sizes.site[[i]][1]-1)])%*%inv.cor.draw[[i]]%*%alpha.draw[k:(k+sizes.site[[i]][1]-1)]
		k <- k+(sizes.site[[i]])[1]
	}
	return(termvalue)
}

invquardtermX <- function(X, d.site, inv.cov.draw){
	sizes.site <- lapply(d.site, dim)
	termvalue <- 0
	k <- 1
	for(i in 1:length(d.site)){
		termvalue <- termvalue+ t(X[k:(k+sizes.site[[i]][1]-1),])%*%inv.cov.draw[(k:(k+sizes.site[[i]][1]-1)), (k:(k+sizes.site[[i]][1]-1))]%*%X[k:(k+sizes.site[[i]][1]-1), ]
		k <- k+(sizes.site[[i]])[1]
	}
	termvalue <- solve(termvalue)
	return(termvalue)
}

invquardtermXY <- function(Y, X, d.site, inv.cov.draw){
	sizes.site <- lapply(d.site, dim)
	termvalue <- 0
	k <- 1
	for(i in 1:length(d.site)){
		termvalue <- termvalue+ t(X[k:(k+sizes.site[[i]][1]-1),])%*%inv.cov.draw[(k:(k+sizes.site[[i]][1]-1)), (k:(k+sizes.site[[i]][1]-1))]%*%Y[k:(k+sizes.site[[i]][1]-1)]
		k <- k+(sizes.site[[i]])[1]
	}
	return(termvalue)
}



detterm <- function(d.site, sizes.noise, inv.cor.draw){
	sizes.site <- lapply(d.site, dim)
	detvalue <- 0
	k <- 1+sizes.noise
	for(i in 1:length(d.site)){
		detvalue <- detvalue+ determinant(inv.cor.draw[[i]])$modulus
		k <- k+(sizes.site[[i]])[1]
	}
	return(detvalue)
}

cov.alpha <- function(d.site, inv.cor.draw, sigma2.alpha.draw, sigma2.draw, sizes.noise, size.points, PID.sizes.list){
	invcov.prior <- vector("list", (length(d.site)+1))
	invcov.prior[[1]] <-diag(1/sigma2.alpha.draw*rep(1,sizes.noise))
	
	invcov.prior[2:(length(d.site)+1)] <- lapply(inv.cor.draw, "/", c(sigma2.alpha.draw) )
	diag.temp <- lapply(mapply(diag, x=PID.sizes.list), "/", sigma2.draw)
	cov.alpha.inv <- mapply("+", invcov.prior, diag.temp)
	cov.alpha <- lapply(cov.alpha.inv, solve)
	return(cov.alpha)
}

mean.alpha.func <- function(d.site,Sigma.alpha,  sigma2.draw, res.fix.group){
	mean.alpha.list <- vector("list", (length(d.site)+1))
	mean.alpha.list[[1]] <- Sigma.alpha[[1]]%*%res.fix.group[1:sizes.noise]/sigma2.draw
	k <- 1+sizes.noise
	for(i in 2:(length(d.site)+1)){
		mean.alpha.list[[i]] <- Sigma.alpha[[i]]%*%res.fix.group[k:(k+sizes.site[[(i-1)]][1]-1)]/sigma2.draw
		k <- k + (sizes.site[[i-1]][1])
	}
	return(mean.alpha.list)		
}

cholcov <- function(Sigma.alpha){
	chol.list <- mapply(chol, Sigma.alpha)
	return(chol.list)
}

alpha.draw.func <- function(alpha.mean, cholcov, sizes.noise, sizes.site){
	alpha.draw.list.sd <- mapply(rnorm, n = append(sizes.noise, lapply(sizes.site, mean)))
	alpha.draw.list <-  mapply("%*%", alpha.draw.list.sd, chol.Sigma.alpha)
	alpha.draw.vec <- t(do.call("cbind", alpha.draw.list)) + (do.call("rbind", alpha.mean))
	return(alpha.draw.vec)
}

incov <- function(d.site, phi.draw, sigma2.alpha.draw, sigma2.draw, sizes.site, size.points){
	d.site.scaled <- lapply(lapply(d.site, "-"), "/", phi.draw) 
	cor <- lapply(d.site.scaled, exp)	
	invcov.stack <- matrix(0, size.points, size.points)
	invcov.stack[(1:sizes.noise),(1:sizes.noise)] <- 1/(sigma2.draw+sigma2.alpha.draw)*diag(rep(1))
	k <- 1
	for(i in 1:length(d.site)){
		invcov.stack[(k:(k+sizes.site[[i]][1]-1)), (k:(k+sizes.site[[i]][1]-1))] <- solve(cor[[i]]*c(sigma2.alpha.draw)+ diag(rep(1, sizes.site[[i]][1]))*c(sigma2.draw))
		k <- k+(sizes.site[[i]])[1]
	}
	return(invcov.stack)
}





krige.incov <- function(d.site,phi.draw, sigma2.alpha.draw, sigma2.draw, inv.sph.cor10.site, size.points){
	d.site.scaled <- lapply(lapply(d.site, "-"), "/", phi.draw) 
	cor <- lapply(d.site.scaled, exp)	
	cov <- lapply(cor, "*", sigma2.alpha.draw)
	sizes.site <- lapply(d.site, dim)
	invcov.stack <- matrix(0, size.points, size.points)
	k <- 1
	for(i in 1:length(d.site)){
		invcov.stack[(k:(k+sizes.site[[i]][1]-1)), (k:(k+sizes.site[[i]][1]-1))] <- solve(cov[[i]]*inv.sph.cor10.site[[i]] + sigma2.draw*diag(rep(1, sizes.site[[i]][1])))
				k <- k+(sizes.site[[i]])[1]
	}
	return(invcov.stack)
}