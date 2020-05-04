#Read the data
x=read.table("../results/x.txt")
y=read.table("../results/y.txt")

#Some parameters
nbodies = 986
nframes = 100

#Generate a vector of colours
colrs = c("#FFFEC5","#FFFBED","#FFFFFF","#E4FEFF","#FFD5DD")


dir.create("../frames")

#Make a big enough blank plot for all the trajectories
for(i in 1:nframes) { #Number of frames, set high to check system evolution quickly
  if(i%%5 == 0) { #Sample every X frames 
    
    png(paste("../frames/frame",i,".png",sep=""),width=500,height=500)
  
    par(bg="black",cex.lab=1.1,col.lab="white",col.axis="white",col.main="white",cex.main=1.8,pch=16)
    plot(x[i,1]-x[i,1],y[i,1]-y[i,1],type="p",cex=1.4,col="slateblue4",xlim=c(-1.5*max(x[1,]),1.5*max(x[1,])),ylim=c(-1.5*max(y[1,]),1.5*max(y[1,])),asp=1,main="Simulación: galaxia.",xlab="x (ndu)",ylab="y (ndu)")
    box(which="plot",lty="solid",col="white")
    for(j in 2:nbodies) { #Plot from the second to the N-th body, centering around the black hole
      points(x[i,j]-x[i,1],y[i,j]-y[i,1],type="p",cex=0.3,col=colrs[j%%5]) #Plots position of the body
    }
    dev.off()
  }
}
