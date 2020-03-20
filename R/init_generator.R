#Please note: it's best to run the whole script to visualize the resulting distribution before
#using the position and velocity data. Some combinations of variables result in a weak
#representation of a spiral galaxy.

N = 15 #Number of spirals
nbodies = 1000 #Approximate number of bodies in the system.
R = 7.72 #Outtermost ellipse radius, or radius of the galaxy
ratio = 0.7 #small semiaxis : big semiaxis ratio

#Normalize R by Sun-Sgt A* distance:
R= R / 2.47

#Generate the ellipses
xaxis = seq(R/N,ratio*R,len=N)
yaxis = seq(R/N,R,len=N)
theta = seq(0,pi/2,len=N) #Angles for 2D rotation
phi = seq(0,2*pi,len=nbodies/N) #Parameter to plot the ellipses

#Initialize position and velocity vectors
x = matrix(NA, ncol = length(phi), nrow = N)
y = vx = vy = x

#Plot the results (only for visual inspection)
plot(0,0,type="n",xlim=c(-R,R),ylim=c(-R,R),asp=1,main="Density wave theory: star distribution.",xlab="x (ndu)",ylab="y (ndu)")
mtext(paste("Number of bodies: ",length(x),sep=""),side=3,line=0.1)

#Generate positions and velocities 
for(i in 1:N) {
  for(j in 1:length(phi)) {
    #Position of each body
    x[i,j] = xaxis[i]*cos(theta[i])*cos(phi[j])-yaxis[i]*sin(theta[i])*sin(phi[j])
    y[i,j] = xaxis[i]*sin(theta[i])*cos(phi[j])+yaxis[i]*cos(theta[i])*sin(phi[j])
    #Velocity of each body
    vx[i,j] = -y[i,j]
    vy[i,j] = x[i,j]
    
    lines(x[i,j],y[i,j],type="p",cex=0.2,col="red")
  }
}

#Normalize the velocities and invert direction
vx = -vx / sqrt(x^2+y^2)
vy = -vy / sqrt(x^2+y^2)

#Now, create a vector of predicted velocities (see notes) and escale vx,vy.
r = sqrt(x^2+y^2)
v=1/r
vx = v*vx
vy = v*vy

#Plot a bunch of velocities to check that nothing weird happened
k = min(dim(x))
for(i in seq(1,k,by=round(k/4))) {
  arrows(x[k,i],y[k,i],x[k,i]+vx[k,i],y[k,i]+vy[k,i],length=0.1,col="black")
  mtext(paste("If you can't see ",k," arrows, some stars are overlapping. Try a smaller N.",sep=""),side=1,line=2,cex=0.8)
  }

#SOME OTHER MAGNITUDES.

#Create a mass vector and normalize it by Sagitarius A's mass.
m = rnorm(length(x),mean=1.0891,sd=0.13)
m = m/(4.143e6)

#Create an effective radius vector (mean = Solar system's radius) and normalize it.
reff = rnorm(N,mean=5.984e12,sd=(5.984e12)/2)
reff = R / (2.469e20)

#Replace first object for central black hole. See notes.
m[1] = 1.0
x[1] = 0.0
y[1] = 0.0
vx[1] = 0.0
vy[1] = 0.0

#Convert each matrix into a single column
x = c(x) + rnorm(length(x),mean=mean(x)/100,sd=sd(x)/1000)
y = c(y) + rnorm(length(y),mean=mean(y)/100,sd=sd(y)/1000)
vx = c(vx) + rnorm(length(vx),mean=mean(vx)/100,sd=sd(vx)/1000)
vy = c(vy) + rnorm(length(vy),mean=mean(vy)/100,sd=sd(vy)/1000)

#Save everything into a data file
fileinfo = data.frame(m,x,y,vx,vy,reff)
write.table(fileinfo,"../dataset.txt",col.names = FALSE, row.names = FALSE, sep = "\t", dec=".")