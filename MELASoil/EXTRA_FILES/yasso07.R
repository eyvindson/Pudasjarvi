yasso07 <-
function(a, t, cl, init, inf, s, z) {
    
    pa = .Fortran("yasso07", a=as.double(a), t=as.double(t),
cl=as.double(cl), init=as.double(init), inf=as.double(inf), s=as.double(s), z=as.double(z) ) 
   
    return(pa$z)
   
}
