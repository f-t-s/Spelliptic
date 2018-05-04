#This module defines the Spelliptic module and the corresponding operator 
#overloading. 
module Spelliptic
using IterativeSolvers
include("SortSparse.jl")
include("Factorization.jl")

export SpelMat
export SpelMatIter
export setCovStat!
export testAccuracy

#A "spelliptic" operator, that will contain the sparse fatorization of an
#elliptic integral operator/covariance matrix.
struct SpelMat
  #Size of the matrix
  N::Int64
  #Number of nonzero entries of the Factor
  nnz::Int64
  #An integer array containing the IDs of the dofs in order of the elimination 
  P::Array{Int64,1}
  #An integer Array that contains as ith entry the number that the dof with Id i
  #has in the elimination ordering
  revP::Array{Int64,1}
  #An array of doubles that contains the "characteristic distance" of each of
  #the dofs, ordered in elimination order
  distances::Array{Float64,1}
  #The upper-triangular matrix containing the sparse factor
  U::SparseMatrixCSC{Float64, Int64}
end

#Constructors
#Constructor to generate a new SpelMat from a d \times N array of N points in 
#d-dimensional space and an accuracy factor \rho.
function SpelMat( x::Array{Float64,2}, rho::Float64 )
  N::Int64 = size(x,2)
  P, revP, distances = sortPoints( x );
  ind, jnd = sparsityPattern( x, P, revP, distances, rho )
  nnz::Int64 = size( jnd, 1 )
  return SpelMat( N,
                  nnz,
                  P,
                  revP,
                  distances,
                  SparseMatrixCSC{Float64, Int64}( N,
                                                   N,
                                                   ind,
                                                   jnd,
                                                   ones(nnz))'
                 )
end

#Constructs a SpelMat, recycling the pattern, by reusing a given P, revP, 
#and rho
function SpelMat( x::Array{Float64,2}, 
                  rho::Float64, 
                  P::Array{Int64,1}, 
                  revP::Array{Int64, 1},
                  distances::Array{Float64, 1} )
  N::Int64 = size(x,2)
  ind, jnd = sparsityPattern( x, P, revP, distances, rho )
  nnz::Int64 = size( jnd, 1 )
  return SpelMat( N,
                  nnz,
                  P,
                  revP,
                  distances,
                  SparseMatrixCSC{Float64, Int64}( N,
                                                   N,
                                                   ind,
                                                   jnd,
                                                   ones(nnz))'
                 )
end

#Operator and function overloading:
#compute the log determinant
import Base.logdet
function logdet( K::SpelMat )
  return 2 * sum( log.( diag( K.U )))
end

#Operator overloading:
import Base.*
function *(K::SpelMat, v )
  return (((K.U * v[K.P])'*K.U)[K.revP])'
end

import Base.*
function *( v, K::SpelMat )
  return (((K.U * v[K.P]')'*K.U)[K.revP])
end

import Base.\
function \( K::SpelMat, v )
  return ((K.U \ ( v[ K.P, : ]' / K.U  )' )[ K.revP, : ] )
end

import Base./
function /( v, K::SpelMat )
  return (( K.U \ ( v[ : , K.P ] / K.U )' )[ : , K.revP ] )'
end

import Base.getindex
function getindex( K::SpelMat, i::Int64, j::Int64 )
  return dot( K.U[ :, K.revP[i]], K.U[ :, K.revP[j]])
end

function getindex( K::SpelMat, ivec, jvec )
  return K.U[ :, K.revP[ivec]]' * K.U[ :, K.revP[jvec]]
end


#Function to set the content of a SpelMat struct to the those given by a 
#stationary covariance function. 
#The variable sigma allows to add an optional nugget term to the covariance.
#Note that for larger nugget terms the approximation property will detoriate,
#Use SpelMatIter when working with large nugget terms.
function setCovStat!( K::SpelMat, x::Array{Float64,2}, covFunc, sigma::Float64 = 0. )
  for i = 1 : K.N
    for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
      K.U.nzval[ j ] = norm( x[ :, K.P[ i ] ] - x[ :, K.P[ K.U.rowval[ j ] ] ] )
    end
  end
  K.U.nzval .= map( covFunc, K.U.nzval )
  if sigma != 0.
    for i = 1 : K.N
      for j = K.U.colptr[ i ] : ( K.U.colptr[ i + 1 ] - 1 )
        if K.U.rowval[ j ] == i
          K.U.nzval[ j ] += sigma 
        end
      end
    end
  end
  @time icholU_high_level!( K.U )
end

function testAccuracy( K::SpelMat, 
                      x::Array{Float64,2}, 
                      covFunc, 
                      NTest::Int64, 
                      sigma::Float64 = 0. )

  testI = randperm( K.N )[ 1 : NTest ]
  testJ = randperm( K.N )[ 1 : NTest ]
  errorMat = zeros( NTest, NTest )
  distMat = zeros( NTest, NTest )
  trueMat = zeros( NTest, NTest )
  appMat = zeros( NTest, NTest )

  for i = 1 : NTest
    for j = 1 : NTest
      trueMat[ i, j ] = covFunc( norm( x[ :, testI[ i ] ] -
             x[ :, testJ[ j ] ] ) )
      if ( testI[i] == testJ[j] )
        trueMat[ i , j ] += sigma
      end
      appMat[ i, j ] = K[ testI[ i ], testJ[ j ]]
      distMat[ i, j ] = norm( x[ :, testI[ i ]] - x[ :, testJ[ j ]] )
    end
  end
  errorMat = abs.( trueMat - appMat )
  return trueMat[:], appMat[:], errorMat[:], distMat[:]
end



###############################################################################
#An extension of SpelMat that treats the "nugget" by using iterative methods.
#The compression is only applied to the part that comes from a Green's function,
#then, conjugate gradient is used to invert the sum of this compressed matrix
#and its nugget.
mutable struct SpelMatIter
  K::SpelMat
  sigma::Float64
  tol::Float64
end

#Constructors
#Constructor to generate a new SpelMatIter from a d \times N array of N points in 
#d-dimensional space and an accuracy factor \rho.
function SpelMatIter( x::Array{Float64,2}, rho::Float64 )
  return SpelMatIter( SpelMat( x, rho ), 0., 1e-3 )
end

#Constructs a SpelMatIter, recycling the pattern, by reusing a given P, revP, 
#and rho
function SpelMatIter( x::Array{Float64,2}, 
                  rho::Float64, 
                  P::Array{Int64,1}, 
                  revP::Array{Int64, 1},
                  distances::Array{Float64, 1} )

  return SpelMatIter( SpelMatIter( x::Array{Float64,2}, 
                      rho::Float64, 
                      P::Array{Int64,1}, 
                      revP::Array{Int64, 1},
                      distances::Array{Float64, 1} ),
                      0., 1e-3 )
end

#Operator and function overloading:
function setCovStat!( K::SpelMatIter, x::Array{Float64,2}, covFunc, sigma::Float64, sigmaimp = 0.)
  setCovStat!( K.K, x, covFunc, sigmaimp )
  K.sigma = sigma
end


#Operator overloading:
import Base.*
function *(K::SpelMatIter, v )
  return (((K.K.U * v[K.K.P])'*K.K.U)[K.K.revP])' + K.sigma * v
end

import Base.*
function *( v, K::SpelMatIter )
  return (((K.K.U * v[K.K.P]')'*K.K.U)[K.K.revP]) + K.sigma * v
end

import Base.\
function \( K::SpelMatIter, v )
  return cg( K, v; verbose = false, tol = K.tol  ) 
end

import Base./
function /( v, K::SpelMatIter )
  return cg( K, v'; verbose = false, tol = K.tol )
end

import Base.getindex
function getindex( K::SpelMatIter, i::Int64, j::Int64 )
  return dot( K.K.U[ :, K.K.revP[i]], K.K.U[ :, K.K.revP[j]]) + 
  K.sigma * convert( Float64, i==j )
end

function getindex( K::SpelMatIter, ivec, jvec )
  return K.K.U[ :, K.K.revP[ivec]]' * K.K.U[ :, K.K.revP[jvec]] +
  K.sigma * speye( K.K.N )[ivec, jvec]
end

import Base.eltype
function eltype( K::SpelMatIter )
  return eltype( K.K.U )
end

import Base.size
function size( K::SpelMatIter, d )
  return size( K.K.U, d )
end

import Base.A_mul_B!
function A_mul_B!( y, K::SpelMatIter, x )
  y .= K * x
end



#the end of Module SpellipticRL
end
