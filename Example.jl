using Distributions
using Distances
include("./Spelliptic.jl")
using Spelliptic
include("./CovFuncs.jl")

#setting the random seed
srand(123)

#Setting problem parameters
d = 2
N = 50000
l = 1.0
nu = 0.5
rho = 4.

#defining covariance function
cov = r -> matern( r, l, nu )

#generating point cloud
dUniform = Uniform(-3,3)
x = rand( dUniform, d, N )

println( "Test without nugget, using SpelMat:" )
#Generating a SpelMat structure
@time K = SpelMat( x , rho ) 
#Setting the content of the spelliptic struct according to the covariance 
#matrix and point cloud:
@time setCovStat!( K, x, cov )

#Measuring elementwise approximation error
trueSol, appSol, error, distance = testAccuracy( K, x, cov, 500 ) 
println( "mean error = ", mean( error ), " median error = ", median( error ) )

vec = randn( N )
#Matrix vector product:
K * vec;
#Solution of linear system
K \ vec
#computing the logdeterminant
println( "logdet(K) = ", logdet( K ))

println( "Test with large nugget, using SpelMatIter:")
sigma = 10.
@time KIter = SpelMatIter( x , rho ) 
@time setCovStat!( KIter, x, cov, sigma )


#Matrix vector product:
KIter * vec;
#Solution of linear system
KIter \ vec

