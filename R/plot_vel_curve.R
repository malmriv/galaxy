#Read the data
x=read.table("../results/x.txt")
y=read.table("../results/y.txt")
vx=read.table("../results/vx.txt")
vy=read.table("../results/vy.txt")
ax=read.table("../results/ax.txt")
ay=read.table("../results/ay.txt")
N = length(x)

#Generate the frames of the animation
dir.create("../curve")
for(i in 1:2000) {
  if(i%%60 == 0) {
  #Define the magnitudes, position and velocity
  pos = as.numeric(sqrt(x[i,1:N]^2+y[i,1:N]^2))
  vel = as.numeric(sqrt(vx[i,1:N]^2+vy[i,1:N]^2))
  
  png(paste("../curve/curve",i,".png",sep=""),width=500,height=500)
  #Plot the graph itself
  par(bg="black",cex.lab=1.1,col.lab="white",col.axis="white",col.main="white",cex.main=1.8,pch=16)
  plot(pos,vel,xlim=c(0,4),ylim=c(0,3),cex=0.3,col="yellow",xlab="distancia (ndu)",ylab="velocidad (nvu)",main="Curva de rotación.")
  box(which="plot",lty="solid",col="white",lwd=1)
  dev.off()
  }
}
