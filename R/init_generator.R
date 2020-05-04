#Author: Manuel Almagro Rivas (malmriv@correo.ugr.es)
#Description: generates a data file readable by a Fortran program. The file
#contains information about where the bodies in a galaxy are, their velocities,
#their masses and some other magnitudes. Most powers of ten have been omitted
#because all data is eventually going to be normalised.

# DESCRIBE THE GALAXY. Approximate number of bodies, radius (in kiloparsecs) and
#mean mass of the stars (in sun masses).
nbodies = 1500
R = 25
meanmass = 1.0891

# KINEMATICS.
# 1) Normalise radius by Sagittarius A* - Sun distance. Generate approximately
#4N/pi Cartesian coordinates in a square of side 2R and keep only those that lie
#in the circumscribed circle of radius R (approximately N bodies).

R = R / 8.178
x = runif(as.integer(4*nbodies/pi),min=-R,max=R)
y = runif(as.integer(4*nbodies/pi),min=-R,max=R)
inside = which(sqrt(x^2+y^2) < R)
x = x[inside]
y = y[inside]
print(length(x))

# 2) Generate versors that are normal to the position vectors, and multiply
#them by their expected speed (Newtonian approximation).

vx = y/(x^2+y^2)
vy = -x/(x^2+y^2)

# MASS DISTRIBUTION.
# Create a mass vector and normalize it by Sagitarius A's mass (in solar masses)
m = rnorm(length(x),mean=meanmass,sd=meanmass/5)
m = m/(4.154e6)

# GRAVITY.
# 1) Create an effective radius vector.
reff = rnorm(length(x),mean=0.001,sd=0.001/3)

# PREPARE THE FILE.
# 1) Replace first object for central black hole. See notes.
m[1] = 1.0
x[1] = 0.0
y[1] = 0.0
vx[1] = 0.0
vy[1] = 0.0
reff[1] = 0.3

#2) Save everything into a data file that Fortran can easily read.
fileinfo = data.frame(m,x,y,vx,vy,reff)
write.table(fileinfo,"../dataset.txt", col.names = FALSE, row.names = FALSE,
sep = "\t", dec=".")
