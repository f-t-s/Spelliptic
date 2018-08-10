include("MutHeap.jl")
#function assignMembership(L::SparseMatrixCSC{Tv,Ti},
#                          distances::Array{Tv,1},
#                          l::Tv,
#                          minBin::Ti) where {Tv<:Real, Ti<:Integer}
#  membership::Array{Ti,1} = Array{Ti,1}(L.n)
#  tempDist::Array{Tv,1} = fill(typemax(Tv),L.n)
#
#  for colIter = 1 : L.n
#    for rowIter = L.colptr[colIter] : (L.colptr[colIter + 1] - 1)
#      j = colIter
#      i = L.rowval[rowIter]
#      val = L.nzval[rowIter]
#      if distances[j] >= l * distances[i]
#        if tempDist[i] > val
#          membership[i] = j
#          tempDist[i] = val
#        end
#      end
#    end
#  end
#  membership[one(Ti)] = one(Ti)
#  return membership
#end

#Defining the subarray Types that are going to be used to store the
#elements of Columns of block sparse matrices and the elements 
#of the blocks of the corresponding ordering.

ColSubArray{Tv,Ti} = SubArray{Tv,2,Array{Tv,2},
                          Tuple{Base.Slice{Base.OneTo{Ti}}, 
                                UnitRange{Ti}},true} where{Tv,Ti}

BlockSparseMatrixCSC{Tv,Ti} = SparseMatrixCSC{ColSubArray{Tv,Ti},
                                              Ti} where{Tv,Ti}

PermSubArray{Ti} =  SubArray{Ti,1,Array{Ti,1},
                             Tuple{UnitRange{Ti}},true} where{Ti}



function superNodalOrder(colptr::Array{Ti,1},
                         rowval::Array{Node{Tv,Ti},1},
                         P::Array{Ti,1},
                         revP::Array{Ti,1},
                         distances::Array{Tv,1},
                         l::Tv,
                         h::Tv,
                         minBin::Ti) where {Tv<:Real,Ti<:Integer}
  #Number of nodes:
  N = size(P,1)
  #Define a helper function that turns a given distance into the 
  #associated level. The larger than distances[2] is considered
  #level 1, between distances[2] and h distances[2] is considered 
  #level 2, and so on.
  function _dist_to_level( r::Tv ) 
    if r <= distances[2]
      return ceil( Ti, log(h,r/distances[2]) ) + one(Ti) + one(Ti)
    else
      return one(Ti)
    end
  end
  #Number of levels, given by the level of the last point:
  q = _dist_to_level( distances[end] )
  #Arrays keeping track of for each point of the bin to which it belongs and
  #for each bin of the number of points belonging to it
  member_count::Array{Ti,2} = fill(zero(Ti),N, q)
  member_ship::Array{Ti,2} = Array{Ti,2}(undef,2,N)
  tempDist::Array{Tv,1} = fill(typemax(Tv),N)

  #Assigns entry a membership:
  for colIter = 1 : N
    for rowIter = colptr[colIter] : (colptr[colIter + 1] - 1)
      #The node with index j is a father of the node with index i.
      j = colIter
      i = rowval[rowIter].id 
      val = rowval[rowIter].val
      #If the node with index j is coarse enough to "fit" index i and it's local
      #neighborhood, it is being considered as a possible cluster.
      if distances[j] >= l * distances[i]
        #We check if the new possible cluster j, which has distance of val 
        #to the point i, is closer to i than the best admissible cluster so far.
        if val < tempDist[i]
          #We add the corresponding membership information to the index i of 
          #membership. Point i is now assigned to the cluster represented by 
          #j,  on the level given by k = l * distances[i].
          #We note, that by definition the point j is part of the first k 
          #levels.
          member_ship[1,i] = _dist_to_level(l * distances[i])
          member_ship[2,i] = j
          tempDist[i] = val
        end
      end
    end
  end
  #Next, we count the number of points put into each bin:
  for index = 1 : N
    member_count[member_ship[2,index], member_ship[1,index]] += one(Ti)
  end
  #Now, starting from the finest level, we want to eliminate supernodes with
  #less than minBin nodes by aggregating them into nodes on coarser scales.

  #TODO: For the simple version we implement so far, we don't worry about this
  #We are sorting the elements from end to start
  #P_sup_elem = sortcols( vcat( member_ship, (one(Ti) : N)' ), rev=true )[3,:]
  for index = one(Ti) : N
    #If an element is in a bin that is too small, we move it up until it 
    #can join a large enough cluster
    newk = member_ship[1,index]
    newInd = member_ship[2,index]
    while member_count[newInd, newk ] < minBin
      #if newk > _dist_to_level( l * distances[newInd] )
      if newk > 1
        newk -= one(Ti)
      elseif newInd > 1
        newk = member_ship[1, newInd ]
        newInd = member_ship[2,newInd]
      else
        break
      end
    end
    #decreasing the member_count at the old location by one
    member_count[ member_ship[2,index], member_ship[1,index] ] -= one(Ti)
    #increasing the member count at the new location
    member_count[ newInd, newk ] += one(Ti)
    #updating the newk and newInd
    member_ship[2,index] = newInd
    member_ship[1,index] = newk
  end

  #Sorting the elements by cluster, from large distance to small distance
  NBin = count( member_count .!= zero(Ti) )
  member_count_lin = cumsum( vcat( 0, vec(member_count)[ 
          (LinearIndices(member_count))[findall(member_count.!= 0) ] ]
                                  )).+ one(Ti)  
  P_sup = Array{PermSubArray{Ti},1}(undef,NBin)
  revP_sup = Array{Ti, 2}(undef,N,2)
  #Permuting with P in order to get  the results with respect to the original
  #ordering, not just the ordering obtained from sortSparse
  P_sup_elem = P[sortslices( vcat( member_ship, (one(Ti) : N)' ), dims=2 )[3,:]]
  for k = 1 : NBin
    P_sup[k] = view(P_sup_elem, member_count_lin[k] : (member_count_lin[k+1] - 1 )
                )
    for j = member_count_lin[k] : (member_count_lin[k+1] - 1 )
      revP_sup[P_sup_elem[j],1] = k 
      revP_sup[P_sup_elem[j],2] = j - member_count_lin[k] + 1 
    end
  end

  return P_sup, revP_sup, P_sup_elem, member_count_lin
end

function superNodalPattern( P_sup::Array{PermSubArray{Ti},1},
                            revP_sup::Array{Ti,2}, 
                            member_count_lin::Array{Ti,1},
                            colptr::Array{Ti,1},
                            rowval::Array{Ti,1},
                            P::Array{Ti,1} ) where {Ti}
N_sup::Ti = size( P_sup, 1 )
N::Ti = size( colptr, 1 ) - one(Ti)


#TODO Can be made much more memory efficient by writing for loops. 
L::SparseMatrixCSC{Bool,Ti} = SparseMatrixCSC{Bool,Ti}( N, N, colptr, rowval,
                             fill(true, size( rowval, 1) ) )
  

  #TODO findn is deprecated in 1.0, to be fixed.
  I::Array{Ti,1} = getindex.(getfield.(findall(!iszero, L), :I), 1)
  J::Array{Ti,1} = getindex.(getfield.(findall(!iszero, L), :I), 2)  
  IJ = unique(hcat( revP_sup[P[I],1], revP_sup[P[J], 1] ),dims=1)
  L = sparse([true])
  I = [zero(Ti)]
  J = [zero(Ti)]
  I = IJ[:,1]
  J = IJ[:,2]
  U_sup::SparseMatrixCSC{Bool,Ti} = sparse( I, J, fill(true, size( I, 1 ) ) )
  I = [zero(Ti)]
  J = [zero(Ti)]
  U_sup = triu( U_sup .| U_sup' )

  return U_sup
end

function superNodalMatrix(P::Array{PermSubArray{Ti},1}, 
                          revP::Array{Ti,2}, 
                          colptr::Array{Ti,1}, 
                          rowval::Array{Ti,1},
                          xfill::Tv) where {Ti,Tv}

  #We have now constructed the sparsity pattern. As the next step, we want to
  #construct the sparser matrix with view entries.
  N = ( size(colptr, 1 ) - 1 )
  nzval =  Array{ColSubArray{Tv,Ti},1}( undef, size(rowval,1) )

  matArray = Array{Array{Tv,2}, 1}(undef, N)
  
  for j = 1 : N
    #We first determine the overall array
    length = zero(Ti)
    for i = colptr[j] : ( colptr[j + 1] - 1 )
      length += size( P[ rowval[ i ] ], 1 )
    end

    #We want the shorter dimension, which is shared by all subarrays in the 
    #column, to come first so that we can have fast indexing for the 
    #subarrays
    matArray[j] = fill( xfill, size( P[j], 1 ), length )
    size( matArray[j] )

    #Now we assign the suitable subarrays to the nonzeroes of 
    startInd::Ti = one(Ti)
    for i = colptr[j] : ( colptr[j + 1] - 1 )
      nzval[i] = 
      view( matArray[j], :, startInd : (startInd + size( P[rowval[ i ] ], 1 ) -
                                        one(Ti)))
      startInd += size( P[ rowval[ i ] ], 1 )
    end
  end

  return BlockSparseMatrixCSC{Tv,Ti}(N,
                                     N,
                                     colptr,
                                     rowval,
                                     nzval)
                        

  
end

function fillMatrix!(U::BlockSparseMatrixCSC{Tv,Ti}, 
                     P_sup::Array{PermSubArray{Ti},1},
                     entryFunc) where{Tv,Ti}
for j_block = 1 : U.n
    for ind_block = U.colptr[j_block] : ( U.colptr[j_block + 1] - 1) 
      i_block = U.rowval[ind_block]
      
      if i_block != j_block
        #The Subarrays representing the individual entries (blocks) are 
        #transposed with respect to the surrounding block matrix.
        #This is to store the external matrix in CSC format while also 
        #allowing for fast indexing of the array views.
        for j = 1 : size(P_sup[i_block], 1)
          for i = 1 : size(P_sup[j_block], 1)
            U.nzval[ind_block][i,j] = entryFunc(P_sup[j_block][i], 
                                              P_sup[i_block][j])
          end
        end
      else
        #First set to zero to then only fill the triangular half
        U.nzval[ind_block] .= zero(Tv)
        for j = 1 : size(P_sup[i_block], 1)
          for i = j : size(P_sup[j_block], 1)
            U.nzval[ind_block][i,j] = entryFunc(P_sup[j_block][i], 
                                              P_sup[i_block][j])
          end
        end
      end
    end
  end
end

function fillMatrix_rbf!(U::BlockSparseMatrixCSC{Tv,Ti}, 
                         P_sup::Array{PermSubArray{Ti},1},
                         x::Array{Tv,2},
                         covFunc) where{Tv,Ti}
  N = size(x,2)
  #Recast as static arrays to use fast methods provided by StaticArrays.jl
  #Possibly remove, doesn't seem to yield a lot.
  xM = reshape( reinterpret(SVector{size(x,1),Tv}, vec(x)), N )
  function entryFunc( i::Ti, j::Ti )
    return covFunc( dot(xM[i],xM[i]) + dot(xM[j],xM[j]) - 2 * dot(xM[i],xM[j]) )
  end
  fillMatrix!( U, P_sup, entryFunc )
end

import Base.full
function full( U::BlockSparseMatrixCSC{Tv,Ti} ) where{Tv,Ti<:Integer}
  #Determine the sizes of the different blocks
  blockSizes = Array{Ti,1}( U.n )
  for k = 1 : U.n
    blockSizes[k] = size( U.nzval[U.colptr[k]],1 )
  end
  #Create array of cumulative block sizes
  cumBlockSizes = cumsum( vcat( zero(Ti), blockSizes ) )

  #Create output matrix
  out = zeros( Tv, cumBlockSizes[U.n + one(Ti) ], cumBlockSizes[U.n + one(Ti)] )
  size(out)

  #Set the matrix entries
  for j_block = 1 : U.n
    for ind_block = U.colptr[j_block] : ( U.colptr[j_block + 1] - 1) 
      i_block = U.rowval[ind_block]
      out[cumBlockSizes[i_block] + 1 : cumBlockSizes[i_block + 1],
          cumBlockSizes[j_block] + 1 : cumBlockSizes[j_block + 1]] = 
      (U.nzval[ ind_block])'
      (U.nzval[ ind_block])'
    end
  end
  return out
end

###############################################################################
#computing the supernodal factorization
###############################################################################
function _updateElem_IChol!( Target, u_nzval, u_rowval, v_nzval, v_rowval )
  u_ind = 1
  v_ind = 1
  #Can be written slightly more simple if it is known which of the two vectors
  #is longer. The minus one is important, since we do not want to update the 
  #k,k th entry usig the k th row.
  s = min(min( size(u_rowval,1), size(v_rowval,1) ), 
          max( size(u_rowval,1), size(v_rowval,1) ) - 1)
  while  max(u_ind, v_ind) <= s
    u_rw = u_rowval[u_ind] 
    v_rw = v_rowval[v_ind] 
    if u_rw == v_rowval[ v_ind ]
      BLAS.gemm!('N', 'T', 
                 -1.0, u_nzval[u_ind], v_nzval[v_ind], 
                 1.0, Target)
      u_ind += 1
      v_ind += 1
    elseif u_rw < v_rw
      u_ind += 1
    else
      v_ind += 1
    end
  end
end

function _updateElem_IChol!(Target, nzval, rowval, 
                            u_min, u_max, 
                            v_min, v_max)
  u_ind = u_min
  v_ind = v_min
  while  ( u_ind <= u_max ) && ( v_ind <= v_max )
    @inbounds u_rv = rowval[u_ind] 
    @inbounds v_rv = rowval[v_ind] 
    if u_rv == v_rv
      @inbounds BLAS.gemm!('N', 'T', 
                 -1.0, nzval[v_ind], nzval[u_ind],
                 1.0, Target)
      #elementUpdate!(nzval[v_ind], nzval[u_ind], Target)
      u_ind += 1
      v_ind += 1
    elseif u_rv < v_rv
      u_ind += 1
    else
      v_ind += 1
    end
  end
end



function IChol!( U::BlockSparseMatrixCSC{Tv,Ti} ) where{Tv, Ti<:Integer}
#  a = view(U.nzval, U.colptr[1] :  (U.colptr[1+1] - 1) )
#  b = view(U.rowval, U.colptr[1] : (U.colptr[1+1] - 1) )
#  c = view(U.nzval, U.colptr[1] : (U.colptr[1+1] - 1) )
#  d = view(U.rowval, U.colptr[1] : (U.colptr[1+1] - 1) )

  for j = 1 : U.n
#    u_nz = view(U.nzval, U.colptr[j] :  (U.colptr[j+1] - 1) )
#    u_rv = view(U.rowval, U.colptr[j] : (U.colptr[j+1] - 1) )
    for i_ind = U.colptr[j] : (U.colptr[j + 1] - 1)
      i = U.rowval[i_ind]
#       v_nz = view(U.nzval, U.colptr[i] : (U.colptr[i+1] - 1) )
#       v_rv = view(U.rowval, U.colptr[i] : (U.colptr[i+1] - 1) )

       _updateElem_IChol!(U.nzval[i_ind], U.nzval, U.rowval,  
                          U.colptr[i], U.colptr[i+1] - 2,
                          U.colptr[j], U.colptr[j+1] - 2)


#@time       _updateElem_IChol!(U.nzval[i_ind], u_nz, u_rv, v_nz, v_rv )


      #Updating the element at question
#      _updateElem_IChol!(U.nzval[i_ind], 
#                         view(U.nzval, U.colptr[j] :  (U.colptr[j+1] - 1) ),
#                         view(U.rowval, U.colptr[j] : (U.colptr[j+1] - 1) ),
#                         view(U.nzval, U.colptr[i] : (U.colptr[i+1] - 1) ),
#                         view(U.rowval, U.colptr[i] : (U.colptr[i+1] - 1) ))

                         
      #Treating the off-diagonal case, by normalising the nonzero with the 
      #corresponding lower triangular part.
      if i < j
      BLAS.trsm!('R', 'L', 'N', 'N', 1.0,
                   U.nzval[ ( U.colptr[ U.rowval[ i_ind ] + 1 ] - 1 ) ],
                   U.nzval[i_ind])
      #Treating the diagonal case:
      else 
        #removing nonzeros on the off-diagonal 
        tril!(U.nzval[ i_ind ])
        #Computing the Cholesky factoriztaion in place, using only the 
        #upper triangular part
        LAPACK.potrf!('L', U.nzval[i_ind])
      end
    end
  end
end

function elementUpdate!(A::ColSubArray{Tv,Ti}, 
                        B::ColSubArray{Tv,Ti}, 
                        C::ColSubArray{Tv,Ti}) where{Tv,Ti<:Integer}
  for j = 1 : size(C,2)
    for k = 1 : size(A,2)
      for i = 1 : size(C,1)
        C[i,j] -= A[i,k] * B[j,k] 
      end
    end
  end
end




