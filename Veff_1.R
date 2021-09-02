###################################################################
# Inputs 
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = 100000
}

r <- 0.35             # stillinger cluster criteria
GRO <- args[1]        # GRO file
L <- args[2]          # Number of MC steps
Np <- rep(NA,L)       # Record Number of particles in each MC step 

###################################################################
# Read gro file into X
library(bio3d)
system("rm -rf PDB.pdb", ignore.stdout = TRUE, ignore.stderr = TRUE)
x <- paste0("echo System | gmx trjconv -f ",GRO, " -s ",GRO, " -o PDB.pdb")
system(x, ignore.stdout = TRUE, ignore.stderr = TRUE)
PDB <- read.pdb("PDB.pdb")
system("rm -rf PDB.pdb", ignore.stdout = TRUE, ignore.stderr = TRUE)

N <- length(PDB$atom$x)     # Number of paricles in gro

X <- matrix(rep(NA,3*(N+1)),nrow = N+1,ncol = 3)
for(i in 1:N)
{
  X[i,] <- cbind(PDB$atom[i,][9:11]$x/10
                 ,PDB$atom[i,][9:11]$y/10
                 ,PDB$atom[i,][9:11]$z/10)
}

###################################################################
# Function to generate random point within a sphere
rsphere <- function(n, r = 0.35, center=cbind(0,0,0)) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  cbind(x+center[,1], y+center[,2], z+center[,3])
}
# Function to calc distance
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

###################################################################
# GCMC 

for(i in 1:L)
{ 
  Np[i] <- (sum(is.finite(X[,1])))         # Count number of particles
  S <- sample(1:N,size = 1,replace = TRUE) # Pick reference particle
  if(Np[i]==N)                             # Always insert if Np=N
  {
  X[(N+1),] <- rsphere(1,center = rbind(X[S,]))
  } else                                   # Try to delete if Np=N+1
  {
    if(euc.dist(X[N+1,],X[S,]) < r)        
    { X[(N+1),] <- cbind(NA,NA,NA)}        # Delete particle if  
  }                                        # within r of ref particle
} 

# 
dG <- -log((4/3)*pi*r^3) - log(sum(Np==(N+1))/sum(Np==(N)))
Veff <- exp( -dG )
cat(Veff)
