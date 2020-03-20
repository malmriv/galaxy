#Read the data
x=read.table("../results/x.txt")
y=read.table("../results/y.txt")


dir.create("../frames")

#Make a big enough blank plot for all the trajectories
for(i in 1:2000) { #Number of frames, set high to check system evolution quickly
  if(i%%4 == 0) { #Sample every X frames in case there are too many of them
    png(paste("../frames/frame",i,".png",sep=""),width=500,height=500)
    par(bg="black",cex.lab=1.1,col.lab="white",col.axis="white",col.main="white",cex.main=1.8,pch=16)
    plot(x[i,1]-x[i,1],y[i,1]-y[i,1],type="p",cex=2,col="slateblue4",xlim=c(-2.0*max(x[1,]),2.0*max(x[1,])),ylim=c(-2.0*max(y[1,]),2.0*max(y[1,])),asp=1,main="Simulación: galaxia.",xlab="x (ndu)",ylab="y (ndu)")
    box(which="plot",lty="solid",col="white")
    for(j in 2:1020) { #Plot from the second to the N-th body, centering around the black hole
      points(x[i,j]-x[i,1],y[i,j]-y[i,1],type="p",cex=0.3,col="khaki1") #Plots position of the planet
    }
    dev.off()
  }
}
