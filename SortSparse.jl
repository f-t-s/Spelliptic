include("./MutHeap.jl")
#############################################################################
#Introducing the "daycare", which keeps track of the descendants of every node
#The struct is essentially a buffered lower triangular CSC sparse matrix,
#together with an of the degrees of freedom
#############################################################################
mutable struct Daycare{Tv,Ti<:Integer}
  NParents::Ti
  NChildren::Ti
  NBuffer::Ti

  #This array gives for contains the ordering. The i-th parent in the 
  #daycare has id P[i]
  P::Array{Ti,1}
  #This array contains as the i-th element the number that the ith parent 
  #has with respect to the multiresolution ordering.
  revP::Array{Ti,1}
  
  #The array that contains the first "child" for every parent
  colptr::Array{Ti,1}

  #The array that contains the global id-s of the children 
  rowval::Array{Node{Tv,Ti},1}
end

#Function that begins a new parent aka column in daycare
function newParent( dc::Daycare{Tv,Ti}, IdParent::Ti ) where 
          {Ti <: Integer,Tv <: Real}
  dc.NParents += 1 
  dc.P[ dc.NParents ] = IdParent
  dc.colptr[ dc.NParents ] = dc.NChildren + 1
  dc.colptr[ dc.NParents + 1 ] = dc.NChildren + 1
  dc.revP[ IdParent ] = dc.NParents
end

function newChild( dc::Daycare{Tv,Ti}, newChild::Node{Tv,Ti}) where 
          {Ti <: Integer, Tv  <: Real}
  if dc.NChildren >= dc.NBuffer
    dc.NBuffer = 2 * dc.NBuffer 
    resize!( dc.rowval, dc.NBuffer )
  end
  dc.NChildren += 1
  dc.colptr[ dc.NParents + 1 ] += 1
  dc.rowval[ dc.NChildren ] = newChild 
end

function newChildren( dc::Daycare{Tv,Ti}, 
                      newChildren::SubArray{Node{Tv,Ti}}) where 
                        {Ti <: Integer, Tv <: Real}
  while dc.NChildren + size(newChildren,1) >= dc.NBuffer - 1
    dc.NBuffer = 2 * dc.NBuffer 
    resize!( dc.rowval, dc.NBuffer )
  end
  dc.NChildren += size(newChildren,1)
  dc.colptr[dc.NParents + 1] += size(newChildren,1)
  dc.rowval[dc.NChildren - size(newChildren,1) + 1 : dc.NChildren] .= newChildren
end

function _determineChildren!(h::MutHeap{Tv,Ti},
                             dc::Daycare{Tv,Ti},
                             parents::Array{Node{Tv,Ti},1},
                             pivot::Node{Tv,Ti},
                             buffer::Array{Node{Tv,Ti},1},
                             rho::Tv,
                             dist2Func) where {Tv<:Real,Ti<:Integer}
#Function to determine the children of a node in the ordering and sparsity
#pattern.
#TODO Update description
#Inputs:
# h:
#   Heap to keep track of the elements and their distance to the points added
#   already.
# dc: 
#   The "daycare" keeping track of the nodes and their children.
# parents:
#   Array that in its i-th position contains the a node with the id of the 
#   preferred parent of the i-th node and the distance of the i-th node to 
#   its preferred parent.
# Id: 
#   The Id of the point, the children of which are being determined.
# rho:
#   The oversampling parameter determining the size of the sparsity pattern.
# dist:
#   The function dist(i,j) gives the distance between the points with ids 
#   i and j.

  #adding the new parent
  newParent( dc, pivot.id )
  dist2ToParent::Tv = parents[pivot.id].val
  lenthscale = pivot.val
  iterBuffer::Ti = zero(Ti)
  for index = dc.colptr[dc.revP[parents[pivot.id].id]] : (dc.colptr[
                    dc.revP[parents[pivot.id].id] + one(Ti)] - one(Ti))
  
    #The candidate point for the pivots children
    candidate::Node{Tv,Ti} = dc.children[ index ]
    #Distance of the candidate to the pivot:
    dist2::Tv = dist2Func( candidate.id, pivot.id )
    #Check whether the candidate is close enough to be added as a child
    if (revP[candidate.id] == zero(Ti)) dist2 <= (lengscale * rho)^2
      dist = sqrt(dist2)
      #Step 1: We add a new child to the pivot:
      #Increase the count of added children by one
      iterBuffer += 1
      #Add node representing the new child to buffer
      buffer[iterBuffer] = Node{Tv,Ti}(dist, candidate.id )
      #Step 2: We possibly update the parent update the distance of the point 
      newDist = update!( h, candidate.id, dist )
      #Step 3: If possible and useful, update the preferred parent:
      if (dist + rho * newDist <= rho * lengscale) && 
         (dist < parents[candidate.id].val)  
         parents[candidate.id] = Node{Tv,Ti}(dist,pivotId)
      end
    #If the candidate is far enough away from the pivots parent, such that it can
    #not possibly be within reach, break:
    elseif (candidate.val - lengthscale * rho)^2 > dist2ToParent
      break
    end
  end
  viewBuffer = view(buffer, 1 : iterBuffer)
  sort!(viewBuffer)
  newChildren(dc, viewBuffer)
end

function sortSparse( N::Ti, rho::Tv, dist2Func, initInd = one(Ti) ) where 
                    {Ti<:Integer, Tv<:Real}
  #Constructing the heap and initialising all variables have maximal distance
  h::MutHeap{Tv,Ti} = MutHeap{Tv,Ti}( Array{Node{Tv,Ti},1}(N), 
                                      one(Ti) : N )
  for i = one(Ti) : N
    h.nodes[i] = Node(typemax(Tv), i)
  end
  #Constructing the Daycare The permutation matrices are initialized to the 
  #maximal element to force errors in case an index is being used before it
  #has been initialized.
  dc::Daycare{Tv,Ti} = Daycare{Tv,Ti}(zero(Ti), 
                        zero(Ti),
                        N,
                        zeros(Ti,N),
                        zeros(Ti,N),
                        zeros(Ti, N + one(Ti)),
                        fill(Node{Tv,Ti}(zero(Tv),zero(Ti)), N))

  #Initializing the Buffer used by _determineChildren!
  nodeBuffer::Array{Node{Tv,Ti},1} = Array{Node{Tv,Ti},1}(N)

  #Initializing the array that will hold the distances measured in the end
  distances::Array{Tv,1} = - ones( Tv, N )

  #Performing the first step of the algorithm:
  #Adding the note initInd as the first parent and making all other nodes its
  #children, as well as updating their distance:
  newParent(dc, initInd)
  distances[1] = typemax(Tv)
  for i = one(Ti) : N
    #adds a new Child and updates the corresponding distance of the child
    newChild( dc, Node{Tv,Ti}(update!(h,i,sqrt(dist2Func(i,initInd))),i))
  end
  parents::Array{Node{Tv,Ti},1} = Array{Node{Tv,Ti},1}(N)
  for i = one(Ti) : N
    parents[i] = Node{Tv,Ti}(sqrt(dist2Func(initInd,i)),i)
  end

  #TODO finish debugging below:
  for i = (one(Ti) + one(Ti) ) : N 
    @show distances[i] = topNode(h).val
    _determineChildren!(h,dc,parents,topNode(h),nodeBuffer,rho,dist2Func)
  end

  return dc.P, dc.revP, distances
end

function sortSparse( x::Array{Tv,2}, rho::Tv, initInd = one(Ti) ) where Tv
  function dist2Func( i::Int64, j::Int64 )
    res::Tv = zero(Tv)
    for d = 1 : size(x,1)
      res += (x[d,i] - x[d,j])^2
    end
    return res
  end
  return sortSparse( size(x,2), rho, dist2Func, initInd )
end



##############################################################################
##Implementation of the ordering algorithm
##############################################################################
#function sortPoints( x::Array{Float64, 2} )
##Function to sort the points in a multicscale way
##Input:
## x:
##   A d times N array containing the measurement points
#  Nbuffer0::Int64 = 10
#  N::Int64 = size( x, 2 )
#  parents::Array{Int64, 1} = Array{Int64, 1}(N)
#  #Create heap:
#  h::Heap = Heap( N, Array{HeapNode, 1}( N ), 1 : 1 : N )
#  for i = 1 : N 
#    h.nodes[ i ] = HeapNode( i, 1000000000000 )
#  end
#  #Create daycare
#  dc::Daycare = Daycare( 0, 0, Nbuffer0, Array{Int64, 1}(N), Array{Int64, 1}(N+1),  Array{Int64, 1}(Nbuffer0), Array{Int64, 1}(N)) 
#  #TODO remove ( only debugging )
#  dc.IdParents .= 0
#  dc.firstChildren .= 0
#  dc.children .= 0
#  dc.lookup .= 0
#  parents .= 0
#  #distances will contain as an i-th entry the distance upon elimnation of the i-th
#  #degree of freedom
#  distances::Array{Float64} = zeros( N ) 
#
#
#  #Setting the first element:
#  #rootID::Int64 = randperm( N )[ 1 ]
#  rootID::Int64 = 1;
#  rootPoint::Array{Float64, 1} = x[ :, rootID ]
#  parents .= rootID
#  newParent( dc, rootID )
#  for i = 1 : N
#    newson( dc, i )
#    distances[1] = max( distances[1], sqrt.( dist2eval( x, i, rootPoint, size( x, 1 ) ) ) )
#    update!( h, i, dist2eval( x, i, rootPoint, size( x, 1 ) ) )
#  end
#  update!( h, rootID, -1. )
#
#  #updating the elements
#  for iter = 2 : N
#    pivotId::Int64 = topId( h )
#    distances[ iter ] = sqrt( topDist2( h ) )
#    _determineChildren!( h, dc, parents, x, pivotId )
#  end
#
#  #Output "P, revP, levels"
#  #Adding a vector with the distances, at which the 
#  return dc.IdParents, dc.lookup, distances
#end
#
##############################################################################
##Implementation of the sparsity pattern
##############################################################################
#function sparsityPattern( x::Array{Float64, 2}, P::Array{Int64, 1}, revP::Array{Int64, 1}, distances::Array{Float64, 1}, rho::Float64 )
#  d::Int64 = size( x, 1 )
#  N::Int64 = size( x, 2 )
#  Nbuffer0::Int64 = 10
#
#  parents::Array{Int64, 1} = Array{Int64, 1}(N)
#  dc::Daycare = Daycare( 0, 0, Nbuffer0, Array{Int64, 1}(N), Array{Int64, 1}(N+1),  Array{Int64, 1}(Nbuffer0), Array{Int64, 1}(N)) 
#
#  dc.IdParents .= 0
#  dc.firstChildren .= 0
#  dc.children .= 0
#  dc.lookup .= 0
#  parents .= 0
#
#  parents .= P[ 1 ]
#  newParent( dc, P[ 1 ] )
#  for i = 1 : N
#    newson( dc, i )
#  end
#
#  for k = 2 : N
#    _determineChildren!( revP, distances, dc, parents, x, P[ k ], rho)
#  end
#  
#  #Sort the jind according to their ordering:
#  for i = 1 : N
#    dc.children[ dc.firstChildren[ i ] : ( dc.firstChildren[ i + 1 ] - 1 ) ] = 
#      P[ sort!( revP[ dc.children[ dc.firstChildren[ i ] : ( dc.firstChildren[ i + 1 ] - 1 ) ] ] ) ]
#  end
#
#  return dc.firstChildren, revP[ dc.children[ 1 : dc.NChildren ] ]
#end
#
