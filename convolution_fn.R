my.convolve <- function(x,y) {
    n <- length(x)
    ny <- length(y)
    Real <- is.numeric(x) && is.numeric(y)
    n1 <- ny - 1
    x <- c(rep.int(0, n1), x)
    n2 <- length(y <- c(rev(y), rep.int(0, n - 1)))
    
    x <- fft(fft(x) * Conj(fft(y)), inverse = TRUE)
    if (Real) return(Re(x)/n2)
    else return(x/n2)
}

my.deconvolve <- function(x,y) {
    n <- length(x)
    ny <- length(y)
    Real <- is.numeric(x) && is.numeric(y)
    n1 <- ny - 1
    x <- c(rep.int(0, n1), x)
    n2 <- length(y <- c(rev(y), rep.int(0, n - 1)))
    
    x <- fft(fft(x) / Conj(fft(y)), inverse = TRUE)
    if (Real) return(Re(x)/n2)
    else return(x/n2)
}

domain <- seq(0,20,0.1)
timepoint <- 1:length(domain)
dist.in1 <- dgamma(domain,2)
dist.in2 <- dexp(domain)
dist.out <- dgamma(domain,3)
convolution <- my.convolve(dist.in1,dist.in2)
convolution.2 <- d(Gammad(2)+Gammad(1))(domain)
scaling <- max(dist.out)/max(convolution)
plot(convolution*scaling)
lines(dist.out)
lines(convolution.2,col="red")

## Try deconvolution
deconvolution <- my.deconvolve(dist.out,dist.in2)
scaling <- max(dist.in1)/max(deconvolution)
plot(deconvolution*scaling)
lines(dist.in1)

## Try it with the cases
plot(my.deconvolve(smooth(cases.symp/sum(cases.symp)),d.inf2test(0:dmax)))
