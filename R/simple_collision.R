#Define the properties of the system
r = 2
N = 17
theta = seq(0,2*pi,len=N+1)
x = r*cos(theta)
y = r*sin(theta)
x = x[-length(x)]
y = y[-length(y)]
m = seq(1,1,length=N)
radius = seq(0.2,0.2,len=N)
vx = rnorm(N,mean=0,sd=0.6)
vy = rnorm(N,mean=0,sd=0.6)


#Arrange the data in a data frame and save into a file
datos = data.frame(m=c(m),x=c(x),y=c(y),vx=c(vx),vy=c(vy),R=c(radius))
write.table(datos,file="../dataset.txt",row.names=FALSE,col.names=FALSE)

#Show the initial setting
plot(x[1],y[1],asp=1,main="Masas.",col=rainbow(N)[1],xlab="x",pch=19,ylab="y",xlim=c(-30,30),ylim=c(-3,3))
for(i in 2:length(x)) {
  lines(x[i],y[i],col=rainbow(N)[i],type="p",pch=19)
}

#EXECUTE THE PROGRAM NOW.

#Read the data
x=read.table("../results/x.txt")
y=read.table("../results/y.txt")

#Check for energy conservation.
energy=read.table("../results/energy.txt")

#Create a folder
dir.create("../frames")
  
#Make a big enough blank plot for all the trajectories
  for(i in 1:dim(x)[1]) { #Number of frames, set high to check system evolution quickly
    if(i%%25 == 0) { #Sample every X frames in case there are too many of them
      png(paste("../frames/frame",i,".png",sep=""),width=500,height=500)
     
       par(bg="black",cex.lab=1.0,col.lab="white",col.axis="white",col.main="white",cex.main=1.8,pch=2.0)
       plot(x[i,1],y[i,1],asp=1,main="N-body simulation.",col=rainbow(N)[1],cex=4.3,pch=19,xlab="x",ylab="y",xlim=c(-3,3),ylim=c(-3,3))
       mtext("Gravity + elastic collisions.",side=3,line=0.1,col="white")
       box(which="plot",lty="solid",col="white")
    
        for(j in 2:N) {
        lines(x[i,j],y[i,j],col=rainbow(N)[j],type="p",pch=19,cex=4.3)
        }
      dev.off()
    }
  }
