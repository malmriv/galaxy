#Import the data:
x = read.table("../results/x.txt")
y = read.table("../results/y.txt")
mass = read.table("../dataset.txt")$V1

#Center around black hole in case there has been any displacement
x = x-x[,1]
y = y-y[,1]

#Define parameters: subdivisions of domain, number of bodies, number of iterations
#and number of iterations to actually plot
subdiv = 20
nbodies = 1517
iterations = 200
analyse = 40

#Define radius to study and actual position of bodies
r_discrete = seq(0,6,len=subdiv)
r_data = sqrt(x^2+y^2)

#Define matrix to save results
density = matrix(ncol=iterations,nrow=subdiv,data=0)

#Read bodies in every ring (every entry of the r_discrete vector) and add their
#mass to the density matrix. Proceed only if there are bodies in the ring (length > 1)
#and avoid considering central black hole bc its mass is too big.
for(j in 1:analyse) {
  index = seq(0,0,len=2)
  for(k in 1:subdiv) {
    index = which(r_data[j,]>=r_discrete[k] & r_data[j,]<r_discrete[k+1])
    if(length(index) > 1) {
      index = sort(index,decreasing=FALSE)
      if(index[1] == 1) index = index[-1]
      density[k,j] = sum(mass[index])
    }
  }
}

image(seq(1,analyse,len=analyse),seq(1,6,len=subdiv),t(density[1:subdiv-1,1:analyse-1]),
col=hcl.colors(100,palette="viridis"),xlab="iteration (adim.)",ylab="radius (ndu)",useRaster=T,main="Radial distribution of mass over time")
