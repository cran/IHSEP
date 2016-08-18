mloglik1b <-
function(jtms,TT=max(jtms),nu,g,Ig=function(x)sapply(x,function(y)integrate(g,0,y)$value),
         Inu=function(x)sapply(x,function(y)integrate(nu,0,y)$value)){
  .External("mloglik1b",jtms,length(jtms),TT,nu,g,Ig,Inu,new.env(),PACKAGE="IHSEP")
}
