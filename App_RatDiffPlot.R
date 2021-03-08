#----------------------#
#-- M&M PLOT for App --#
#----------------------#

# the following is the function for generating the ratio M&M plot
RatPlot <- function(xlim=c(0,20),ylim=c(0,20), # range of the data
                    xlab=expression(paste(m[C]," = median OS for C")),
                    ylab=expression(paste(m[Rx]," = median OS for Rx")), # label of two axes
                    med.surv.Rx=c(10,18,14), # median survival time for the Rx (sub)groups
                    med.surv.C=c(8,10,9), # median survival time for the C (sub)groups
                    # radius=c(sqrt(20),sqrt(25),sqrt(30)), # the length of the arcs
                    angle.1=c(1.65,2.5,2.16), # angle.1 (a vector) = the upper bounds of simultaneous CIs for the ratios
                    angle.2=c(0.85,1.1,0.96), # angle.2 (a vector) = the lower bounds of simultaneous CIs for the ratios
                    color=c(2,3,4), # the color of the lines/arcs
                    pch=c(21,22,23), # the shape of points (to represent the median survival times)
                    legend=c(": g+ group",": g- group", ": {g+,g-} combined group"), # the legend of the plot
                    axis.1=c(expression(paste(m[C],"+")),expression(paste(m[C],"-")), expression(paste(m[C],"*"))),
                    axis.2=c(expression(paste(m[Rx],"+")),expression(paste(m[Rx],"-")), expression(paste(m[Rx],"*")))# labels of the two axes
){
  
  # calculate the median survival ratios for the (sub)groups  			
  ratio = med.surv.Rx/med.surv.C
  
  par(mar=c(5, 5, 2, 2), cex.lab=1.2, cex.axis=1.2)
  # generate a blank plot with axes labels
  plot(xlim,ylim,type="n",xaxt='n',yaxt='n',xaxs="i",yaxs="i",
       ylab=ylab,
       xlab=xlab,asp=1)
  
  # add the 45 degree diagonal line
  abline(a=0,b=1,lwd=1,lty=2,col="grey") 
  # dots for the observed median survival
  lines(med.surv.C,med.surv.Rx,type="p",lwd=1.5,pch=pch,col="black",bg=color)
  for(i in 1:length(ratio)){
    # add lines and arcs
    abline(a=0,b=ratio[i],lwd=1.5,lty=1,col=color[i])
    draw.arc(0,0,radius=sqrt(med.surv.C[i]^2+med.surv.Rx[i]^2),
             angle1=atan(angle.1[i]),angle2=atan(angle.2[i]),n=0.05,col=color[i],
             lwd=1.5,xaxs="i",yaxs="i",lend=1)
    
    # guidelines (vertical/horizontal) for the observed median survival
    lines(c(med.surv.C[i],med.surv.C[i]),c(0,med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
    lines(c(0,med.surv.C[i]),c(med.surv.Rx[i],med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
  }
  
  #axis labels
  axis(1, at=med.surv.C, lab=axis.1)
  axis(2, at=med.surv.Rx, lab=axis.2)
  
  #legend
  legend("topright", legend, cex=1.2,col="black",pt.bg=color, pch=pch)
  
}


DifPlot <- function(xlim=c(0,20),ylim=c(0,20), # range of the data
                    xlab=expression(paste(m[C]," = median OS for C")),
                    ylab=expression(paste(m[Rx]," = median OS for Rx")), # label of two axes
                    med.surv.Rx=c(10,18,14), # median survival time for the Rx (sub)groups
                    med.surv.C=c(8,10,9), # median survival time for the C (sub)groups
                    # radius=c(sqrt(20),sqrt(25),sqrt(30)), # the length of the arcs
                    UCI=c(3,10, 7), # UCI (a vector) = the upper bounds of simultaneous CIs for the differences
                    LCI=c(1, 6, 3), # LCI (a vector) = the lower bounds of simultaneous CIs for the differences
                    color=c(2,3,4), # the color of the lines/CI
                    pch=c(21,22,23), # the shape of points (to represent the median survival times)
                    legend=c(": g+ group",": g- group", ": {g+,g-} combined group"), # the legend of the plot
                    axis.1=c(expression(paste(m[C],"+")),expression(paste(m[C],"-")), expression(paste(m[C],"*"))),
                    axis.2=c(expression(paste(m[Rx],"+")),expression(paste(m[Rx],"-")), expression(paste(m[Rx],"*"))) # labels of the two axes
){
  # calculate the differences of median survival for the (sub)groups				
  dif  = med.surv.Rx-med.surv.C
  
  par(mar=c(5, 5, 2, 1), cex.lab=1.2, cex.axis=1.2)
  # generate a blank plot with axes labels
  plot(xlim,ylim,type="n",xaxt='n',yaxt='n',xaxs="i",yaxs="i",
       ylab=ylab,
       xlab=xlab,asp=1)
  
  # add the 45 degree diagonal line
  abline(a=0,b=1,lwd=1,lty=2,col="grey") 
  # dots for the observed median survival
  lines(med.surv.C, med.surv.Rx,type="p",lwd=1.5,pch=pch,col="black",bg=color)
  
  for(i in 1:length(dif)){
    # add lines for prognostic and predictive
    abline(a=dif[i],b=1,lwd=1.5,lty=1,col=color[i])
    
    # CI for the difference
    lines(c((med.surv.C[i]-(UCI[i]-dif[i])/2), (med.surv.C[i]+(dif[i]-LCI[i])/2)),
          c((med.surv.Rx[i]+(UCI[i]-dif[i])/2), (med.surv.Rx[i]-(dif[i]-LCI[i])/2)), type='l', lwd=1.5, col=color[i], lend=1)
    
    # guidelines for the observed median survival
    lines(c(med.surv.C[i],med.surv.C[i]),c(0,med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
    lines(c(0,med.surv.C[i]),c(med.surv.Rx[i],med.surv.Rx[i]),type="l",lwd=1,lty=3,col=color[i])
  }
  
  #axis labels
  axis(1, at=med.surv.C, lab=axis.1)
  axis(2, at=med.surv.Rx, lab=axis.2)
    
  #legend
  legend("topright", legend, cex=1.2,col="black",pt.bg=color, pch=pch)
  
}