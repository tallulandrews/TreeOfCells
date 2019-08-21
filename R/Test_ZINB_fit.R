
check.ZINB.fit <- function(x, lib.size=rep(1, length(x)), max_r=10^10, e=0.00001) {
	# Fit
	lib.size <- lib.size/mean(lib.size)
	d.obs <- mean(x==0)
	if (d.obs == 0) {return(check.NB.fit(x, lib.size))}
	d_curr <- d.obs;
        d_prev <- -100;
	nc <- length(x);
	while( abs(d_curr-d_prev) > e ) {
		mus <- sum(x)/(nc-d_curr*nc)*lib.size
		weights <- rep(1, times=length(x))
		weights[x == 0] <- (1-d_curr/d.obs)
		obs_err <- sum( (x - mus)^2*weights )
		rg <- sum( mus^2*weights )/(obs_err - sum(mus*weights))
                if (rg <= 0) {rg <- max_r}
		pds <- (1 + mus/rg)^(-rg)
                d.exp <- mean(pds)
		d_prev <- d_curr
                d_curr <- (d.obs - d.exp)
		if (d_curr <= 0) {d_curr <- d_prev}
	}
	# params : mus, d_prev, rg
	p0s <- d_prev+sapply(mus, function(m) {pnbinom(0, mu=m, size=rg)})
	ps <- sapply(1:length(x), function(i) {pnbinom(x[i], mu=mus[i], size=rg)})
	ps[x==0] <- p0s[x==0];
	ll <- log10(prod(ps))
	return(ll);
}

