#TWIG: a function for visualizing the results of a (Julia-based) TWIG analysis
#  eigpath: the location of the eigvec_[abr].txt and eigval_[abr].txt (and tmaxs)
#  abr: the abbreviation of the bifurcation being analyzed, e.g. pitchf
#  pars: an expression listing the parameters and their values, e.g. expression(y[0]==1.2, r==0)
#  npart: the number of participation factor bar graphs to display at the top of the rainbow plot
#  formu: an expression of the bifurcation formula to display on the rainbow plot (not shown in if NA)
#  paramlbl: whether or not to label the important parameters on the participation factors bar graphs
twig<-function(eigpath, abr, pars, npart=2, formu=NA, paramlbl=T) { 
  tocol<-function(ix) {
    vv=as.character(as.hexmode(round(ix*255)))
    vv=paste0(ifelse(nchar(vv)==1,"0",""),vv)
    toupper(vv)
  }
  
  layout(matrix(0:npart+1,nc=1),height=c(rep(.12,npart),1-.12*npart))
  npar=length(pars)
  
  par(mar=c(0,3,0,.5),mgp=c(1.8,.7,0),cex=1.1)
  rv<-read.table(paste0(eigpath,'eigvec_',abr,'.txt')) #eigvec
  ntv=ncol(rv)/npar
  for(j in 1:npart) {
    barplot((pfj=as.matrix(rv[,0:(ntv-1)*npar+j]^2)),space=0,names.arg=rep("",ntv),ylab="Particip.",col=rainbow(npar),axes=F,xlim=c(.5,ntv-.5))
    axis(2,at=c(.2,.5,.8))
    apply(pfj,2,which.max)->ldpar; apply(pfj,2,max)->ldval #find streaks where one parameter has participation >2/3 
    ifelse(ldval>2/3,ldpar,NA)->hasdom
    hasdom[is.na(hasdom)]<- 0
    ix0=which.min(is.na(hasdom)); nstrk=1 #first non-NA
    while(ix0+nstrk<length(hasdom)) {
      if(hasdom[ix0]==hasdom[ix0+nstrk]) nstrk=nstrk+1
      else {
        if(hasdom[ix0]>0) text(ix0+nstrk/2-.5,.5,parse(text=gsub(" ?==.+$","",pars[hasdom[ix0]])),cex=1.5) #label that streak
        ix0=ix0+nstrk; nstrk=1
      }
    }
    if(hasdom[ix0]>0) text(ix0+nstrk/2-.5,.5,parse(text=gsub(" ?==.+$","",pars[hasdom[ix0]])),cex=1.5) #label terminal streak
  }
  
  par(mar=c(3,3,0,.5),mgp=c(1.8,.7,0),cex=1.1)
  z<-read.table(paste0(eigpath,'eigval_',abr,'.txt'))
  if(file.exists(paste0(eigpath,'tmaxs_',abr,'.txt')))
    tmax=read.table(paste0(eigpath,'tmaxs_',abr,'.txt'))[,1]
  else tmax=read.table(paste0(eigpath,'tmaxs.txt'))[,1]
  print(dim(z)); print(length(tmax))
  plot(tmax,unlist(z[1,]),type="n",log='xy',ylim=c(max(1e-16,min(z)),max(z)),xlab=expression(t[max]),
       ylab=expression(lambda),axes=F)
  abline(h=1,v=1); box()
  axis(1,at=10^(-2:3),labels=expression(10^-2,10^-1,1,10,100,1000))
  axis(2,at=10^(-8:1*2),labels=expression(10^-16,10^-14,10^-12,10^-10,10^-8,10^-6,10^-4,10^-2,1,100))
  #grid()
  for(j in 1:npar) {
    lines(tmax,unlist(z[j,]),lwd=2)
    for(k in 1:npar)
      points(tmax,z[j,],col=paste0(rainbow(npar)[k],tocol(unlist(rv[k,0:(ntv-1)*npar+j]^2))),pch=16,cex=2)
  }
  legend("bottom",legend=pars,fill=rainbow(npar),cex=.95,bg="white") #title=ti[i],
  if(!is.na(formu)) text(10^par()$usr[1],10^par()$usr[4]/10,formu,pos=4)
}

#twig("~/bifurcations/","pitchf", pars=expression(y[0]==.5,r==0,alpha[1]==0,alpha[2]==0,alpha[3]==0,alpha[4]==0,alpha[5]==0),
# npart=2, formu=expression(dot(y)==r*y-y^3+alpha[1]*y^4+alpha[2]*y^5+alpha[3]*y^6+alpha[4]*y^7+alpha[5]*y^8))