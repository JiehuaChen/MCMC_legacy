source("main_data_Clay.R")
source("invcor_list_v2.R")
source("likfunc_nugget.R")


library(MASS)


##start of MCMC
nIter <- 5000
nBurnin <- 2500
thin <- 10

#MCMC iterations
Y <- Y/100
XTX <- t(X)%*%X
size.points.cluster <- sum(table.dbscan[-1])
size.points <- sum(table.dbscan)
sizes.site <- lapply(d.site, dim)
cluster.asign.order264ed <- sort(cluster.asign)
PID.sizes <- as.vector(table(PIDn))
PID.sizes.list  <-  vector("list", (length(d.site)+1))
sizes.noise <- table.dbscan[1]

PID.sizes.list[[1]] <- PID.sizes[1:sizes.noise]
k <- sizes.noise+1
for(i in 1:(length(d.site))){
	PID.sizes.list[[(i+1)]] <- PID.sizes[(k:(k+sizes.site[[i]][1]-1))]
	k <- k+(sizes.site[[i]])[1]
}

#posterior results
K <- length(unique(PIDn))

n.chains <- 3
keepn <- ceiling((nIter-nBurnin)/thin)
beta.post <- matrix(NA, keepn, dim(X)[2])
alpha.post <- matrix(NA, keepn, K)
sigma2.post <- rep(NA, keepn)
sigma2.alpha.post <- rep(NA, keepn)
phi.post <- rep(NA, keepn)
lik.post <- rep(NA, keepn)
npara.total <- dim(cbind(beta.post, alpha.post, sigma2.post, sigma2.alpha.post, phi.post))[2]
mcmc.results <- array(NA, c(keepn, n.chains, npara.total))

# priors:
# beta  ~ N(0, 0.00001)
# sigma ~ scaled-inv-chisq(0.0001, 1)
v0 <- 0.0001
sigma02 <- 1

#the prior of phi is uniform between [a, b]
a <- 0.01
b <- taper.range*5

for(j in 1:n.chains){
	#inital values
	lm.XY <- lm(Y~X-1)
	beta_ini <- lm.XY$coefficients

	alpha_ini <- rep(0, N)
	alpha.draw <- alpha_ini
	alpha.draw.full <- alpha.draw[PIDn]
	linkfunc.hat <- X%*%beta_ini + alpha.draw
	sigma2_alpha_ini <- runif(1)
	sigma2.alpha.draw  <- sigma2_alpha_ini
	sigma2_ini <- runif(1)
	sigma2.draw <- sigma2_ini
	repeat{
		linkfunc.draw <- X%*%beta_ini + rnorm(N, 0, sqrt(sigma2.draw))
		logit.tmp <- exp(linkfunc.draw)/(1+exp(linkfunc.draw))
		if(sum(logit.tmp[Y<1]==1)==0){
			break
		}		
	}
	
	linkfunc.jump <- rep(sigma2.draw, N)

	phi.draw <- runif(1, a, b)
	inv.cor.draw  <- incor(d.site, phi.draw, sph.cor20.site)
	#variance of the uniform prior of phi
	phi.jump <- (b-a)^2/12

	#index for keep values
	k <- 1
 	#index for jump counts
 	p.ct <- 0
 	last.50p <- matrix(NA, 50, N+1) 	
	ratio.draw <- 1
	
	for(i in 1:nIter){
		ptm <- proc.time()
		#draw link function
		linkfunc.old <- linkfunc.draw
		ratio.old <- ratio.draw
		logit.old <- exp(linkfunc.draw)/(1+exp(linkfunc.draw))
		shape1 <- logit.old*ratio.old 
		shape2 <- (1-logit.old)*ratio.old
		
		lik.link.old <- log(dbeta(Y, shape1, shape2)) + (dnorm(linkfunc.old, mean = linkfunc.hat, sd =sqrt(sigma2.draw), log = TRUE))
		linkfunc.new <- linkfunc.draw + rnorm(N, 0, sqrt(linkfunc.jump))		
		logit.new <- exp(linkfunc.new)/(1+exp(linkfunc.new))
		shape1 <- logit.new*ratio.draw 
		shape2 <- (1-logit.new)*ratio.draw 
		lik.link.new <- log(dbeta(Y, shape1, shape2)) + (dnorm(linkfunc.new, mean = linkfunc.hat, sd =sqrt(sigma2.draw), log = TRUE))
		prob.lik.diff <- lik.link.new - lik.link.old
		prob.lik.diff <- ifelse(is.nan(prob.lik.diff)|is.infinite(prob.lik.diff), 0, prob.lik.diff)
		jump.link <-  rbinom(N, 1, exp(pmin(prob.lik.diff, 0)))
		linkfunc.draw <- linkfunc.new*jump.link + (1-jump.link)*linkfunc.old
		linkfunc.p <- exp(pmin(prob.lik.diff, 0))

		
		res.randef <- linkfunc.draw - alpha.draw.full
		XTR <- t(X)%*%res.randef
	
		#draw fixed effect
		beta.mean <- solve(XTX, XTR)
		beta.draw <- mvrnorm(1, beta.mean, solve(XTX)*sigma2.draw)

		#draw random effect		
		res.fix <- linkfunc.draw- X%*%as.matrix(beta.draw)	
		res.fix.group <-(aggregate(res.fix, list(PID = PIDn), sum))[,2]							
		Sigma.alpha <- cov.alpha(d.site, inv.cor.draw, sigma2.alpha.draw, sigma2.draw, sizes.noise, size.points, PID.sizes.list)
		chol.Sigma.alpha <- cholcov(Sigma.alpha)
		alpha.mean <- mean.alpha.func(d.site, Sigma.alpha, sigma2.draw, res.fix.group)
		alpha.draw <- alpha.draw.func(alpha.mean, chol.Sigma.alpha, sizes.noise, sizes.site)
		
		#draw phi
		#the prior of phi is uniform between [a, b]	
		phi.old <- phi.draw
		inv.cor.draw.old <- inv.cor.draw
		quardterm.value <- quardterm(d.site, alpha.draw, sizes.noise, inv.cor.draw.old)
		lik.old <- -1/(2*sigma2.alpha.draw)*quardterm.value+ 1/2*detterm(d.site, sizes.noise, inv.cor.draw.old)
			
		phi.new <- rnorm(1, phi.old, sqrt(phi.jump))
		
		if(phi.new>=b||phi.new<=a){
			jump.phi <- 0
			phi.p <- 0
			phi.draw <- phi.old
			inv.cor.draw <- inv.cor.draw.old
		}else{
			inv.cor.draw.new <- incor(d.site, phi.new, sph.cor20.site)
			quardterm.value <- quardterm(d.site, alpha.draw, sizes.noise,inv.cor.draw.new)
			lik.new <- -1/(2*sigma2.alpha.draw)*quardterm.value + 1/2*detterm(d.site, sizes.noise, inv.cor.draw.new)
			prob.lik.diff <- lik.new-lik.old
			jump.phi <-  rbinom(1, 1, exp(pmin(prob.lik.diff, 0)))
			phi.draw <- phi.new*jump.phi + (1-jump.phi)*phi.old
			if(jump.phi==1){
			inv.cor.draw <- inv.cor.draw.new}
			else{inv.cor.draw <- inv.cor.draw.old }
			phi.p <- exp(pmin(prob.lik.diff, 0))
		}
		print(phi.draw)
		print(phi.jump)	

		p.ct <- p.ct+1
		last.50p[p.ct, ] <- c(linkfunc.p, phi.p)
		if(p.ct==50){
			p.ct <- 0
		#	jump.d <- c(beta.jump, alpha.jump, phi.jump,error.jump)
			jump.d <- c(linkfunc.jump, phi.jump)
			p.mean <- colMeans(last.50p)
			jump.d <- pmax(pmin(log(0.4)*jump.d/log(p.mean), 10*jump.d), 0.0001)
			linkfunc.jump <- jump.d[1:N]
			phi.jump <- jump.d[(N+1)]
		}		
		#draw standard deviations
		alpha.draw.full <- alpha.draw[PIDn]	
		res.sqsum <- sum((res.fix-alpha.draw.full)^2)
		scale.draw <- (v0*sigma02 + res.sqsum)
		sigma2.draw <- 1/rchisq(1, df=v0+N)*scale.draw
	
		sqsum.alpha <- sum((alpha.draw[1:sizes.noise])^2) + quardterm(d.site, alpha.draw, sizes.noise,inv.cor.draw)
		scale.draw <- (v0*sigma02 + sqsum.alpha)
		sigma2.alpha.draw <- 1/rchisq(1, df=v0+K)*scale.draw
		
		if(i>nBurnin&&((i-nBurnin)%%thin==0)){
			beta.post[k, ] <- beta.draw
			alpha.post[k, ] <- alpha.draw
			sigma2.post[k] <- sigma2.draw
			sigma2.alpha.post[k] <- sigma2.alpha.draw
			phi.post[k] <- phi.draw
			k <- k+1
			print(k)
		}
		print(i)
		print(proc.time()-ptm)

		}
	mcmc.results[,j,] <- cbind(beta.post, alpha.post, sigma2.post, sigma2.alpha.post, phi.post)
}

save("mcmc_results_Clay.RData")

