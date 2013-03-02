log.likfunc.binom <- function(Y, X, beta.draw, alpha.draw, sigma2.alpha.draw, inv.cor.draw, error.draw){	
	linkfunc <- 	X%*%beta.draw + alpha.draw + error.draw
	log.likfunc.fixed <- (t(Y)%*%linkfunc) - sum(4*log(1 + exp(linkfunc)))
	log.likfunc.random <- -1/(2*sigma2.alpha.draw)*t(alpha.draw)%*%inv.cor.draw%*%alpha.draw+ 1/2*determinant(inv.cor.draw)$modulus -N/2*sigma2.alpha.draw
	log.likfunc <- log.likfunc.fixed + log.likfunc.random
	return(log.likfunc)
}

log.likfunc.binom.linkfunc <- function(Y, linkfunc){
	log.likfunc.fixed <-  ((Y)*linkfunc) - (4*log(1 + exp(linkfunc)))
	return(log.likfunc.fixed)
}

log.likfunc.binom.fixed <- function(Y, X, beta.draw, alpha.draw, error.draw){
	linkfunc <- 	X%*%beta.draw + alpha.draw + error.draw
	log.likfunc.fixed <-  (t(Y)%*%linkfunc) - sum(4*log(1 + exp(linkfunc)))
	return(log.likfunc.fixed)
}

log.likfunc.binom.random <- function(Y, fixeffect, alpha.draw, sigma2.alpha.draw, inv.cor.draw, d.site, error.draw){	
	linkfunc <- fixeffect + alpha.draw + error.draw
	log.likfunc.fixed <- (t(Y)%*%linkfunc) - sum(4*log(1 + exp(linkfunc)))
	quardterm.value <- quardterm(d.site, alpha.draw, inv.cor.draw)
	log.likfunc.random <- -1/(2*sigma2.alpha.draw)*quardterm.value
	log.likfunc <- log.likfunc.fixed + log.likfunc.random
	return(log.likfunc)
}

log.likfunc.binom.random.site <- function(Y, fixeffect, alpha.draw.site, sigma2.alpha.draw, inv.cor.draw.site, d.site, error.draw.site){	
	linkfunc <- 	fixeffect + alpha.draw.site + error.draw.site
	log.likfunc.fixed <- (t(Y)%*%linkfunc) - sum(4*log(1 + exp(linkfunc)))
	log.likfunc.random <- -1/(2*sigma2.alpha.draw)*t(alpha.draw.site)%*%inv.cor.draw.site%*%alpha.draw.site
	log.likfunc <- log.likfunc.fixed + log.likfunc.random
	return(log.likfunc)
}

log.likfunc.binom.error <- function(Y, freffect, error.draw, sigma2.draw){	
	linkfunc <- freffect + error.draw
	log.likfunc.fixed <-  Y*linkfunc - (4*log(1 + exp(linkfunc)))
	log.likfunc.error <- -1/(2*sigma2.draw)*(error.draw^2)
	log.likfunc <- log.likfunc.fixed + log.likfunc.error
	return(log.likfunc)
}

log.likfunc.binom.new <- function(Y, X, beta.draw, alpha.draw, sigma2.alpha.draw, sigma2.draw, inv.cor.draw, d.site, error.draw){	
	linkfunc <- 	X%*%beta.draw + alpha.draw + error.draw
	log.likfunc.fixed <- (t(Y)%*%linkfunc) - sum(4*log(1 + exp(linkfunc)))
	quardterm.value <- quardterm(d.site, alpha.draw, inv.cor.draw)
	log.likfunc.random <- -1/(2*sigma2.alpha.draw)*quardterm.value+ 1/2*detterm(d.site, inv.cor.draw)-N/2*sigma2.alpha.draw-N/2*sigma2.draw-sum(error.draw^2)
	log.likfunc <- log.likfunc.fixed + log.likfunc.random
	return(log.likfunc)
}
