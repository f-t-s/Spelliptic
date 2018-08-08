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
  member_ship::Array{Ti,2} = Array{Ti,2}(2,N)
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
  member_count_lin = cumsum( vcat( 0, vec(member_count)[ find(member_count) ]
                                  )).+ one(Ti)  
  P_sup = Array{SubArray{Ti,1,Array{Ti,1},Tuple{UnitRange{Ti}},true},1}(NBin)
  revP_sup = Array{Ti, 2}(N,2)
  P_sup_elem = P[sortcols( vcat( member_ship, (one(Ti) : N)' ) )[3,:]]
  for k = 1 : NBin
    P_sup[k] = view(P_sup_elem, member_count_lin[k] : (member_count_lin[k+1] - 1 )
                )
    for j = member_count_lin[k] : (member_count_lin[k+1] - 1 )
      revP_sup[j,1] = k 
      revP_sup[j,2] = j - member_count_lin[k] + 1 
    end
  end

  return P_sup, revP_sup, P_sup_elem, member_count_lin
end

function superNodalPattern( P_sup::Array{SubArray{Ti,1,Array{Ti,1},
                                    Tuple{UnitRange{Ti}},true},1},
                            revP_sup::Array{Ti,2}, 
                            member_count_lin::Array{Ti,1},
                            colptr::Array{Ti,1},
                            rowval::Array{Ti,1} ) where {Ti}
N_sup::Ti = size( P_sup, 1 )
N::Ti = size( colptr, 1 ) - one(Ti)


#TODO Can be made much more memory efficientby writing for loops. 
L::SparseMatrixCSC{Bool,Ti} = SparseMatrixCSC{Bool,Ti}( N, N, colptr, rowval,
                             fill(true, size( rowval, 1) ) )
  

  I::Array{Ti,1}, J::Array{Ti,1} = findn( L )   
  #Deallocating rowval
  #rowval = [rowval[1]]
  gc()
  IJ = unique(revP_sup[ hcat( I,J ),1],2)
  L = sparse([true])
  I = [zero(Ti)]
  J = [zero(Ti)]
  gc()

  I = IJ[:,1]
  J = IJ[:,2]
  gc()
  U_sup::SparseMatrixCSC{Bool,Ti} = sparse( I, J, fill(true, size( I, 1 ) ) )
  I = [zero(Ti)]
  J = [zero(Ti)]
  gc()
  U_sup = triu( U_sup .| U_sup' )

  #We have now constructed the sparsity pattern. As the next step, we want to
  #construct the sparser matrix with view entries.
#  nnz =  Array{SubArray{Ti,1,Array{Tv,2},Tuple{Colon,UnitRange{Ti}},true},1}
#
#  matArray = Array{Array{Tv,2}, 1}(N_sup)
#  
#
#  for k = 1 : N_sup
#    matArray[k] = Array{}
#    k
#  end


  return U_sup
end

