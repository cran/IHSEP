simPois <-
    function(int=NULL,cens=1,
             int.M=optimize(int,c(0,cens),maximum=TRUE)$obj*1.1,
             B=max(as.integer(sqrt(int.M*cens)),10),
             exp.int=FALSE,
             par=c(1,2)
             ) {
        if(!exp.int){
            tms <- rexp(as.integer(int.M*cens+2*sqrt(int.M*cens)),rate=int.M)
            while(sum(tms)<cens)tms<-c(tms,rexp(B,rate=int.M))
            cumtms <- cumsum(tms)
            tms <- cumtms[cumtms<=cens]
            N <- length(tms)
            if(N==0)return(numeric(0))
            tms[as.logical(mapply(rbinom,n=1,size=1,
                                  prob=int(tms)/int.M) ) ]
        }else{
            P <- pexp(cens,rate=par[2])
            NT <- rpois(1, lambda=par[1]/par[2] * P)
            U <- runif(n = NT, min = 0, max = 1)
            raw <- -log1p(-P*U)/par[2]
            sort(raw)
        }
    }
