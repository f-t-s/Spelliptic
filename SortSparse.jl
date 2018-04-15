#############################################################################
#This file contains the implementation of the the functions that create the 
#multiresolution ordering and sparsity pattern from an input consisting of 
#a set of points
#############################################################################



#############################################################################
#Implementation of the Heap
#############################################################################

#HeapNode containing the global ID and the distance squared to the points added
mutable struct HeapNode
  Id::Int64
  dist2::Float64
end

#Heap containing its size (for convenience), an array of HeapNodes, 
#and and array mapping the Id's of the nodes to their position in the heap.
struct Heap
  N::Int64
  nodes::Array{HeapNode,1}
  lookup::Array{Int64,1}
end

#Funktion swapping the nodes a and b, as identified by their respective heap ID,
#in the heap. Presently implemented by swapping the two values of the nodes.
function _swap!( h::Heap, aheapInd::Int64, bheapInd::Int64 )
  #Assinging the corresponding new values to the lookup table
  h.lookup[ h.nodes[aheapInd].Id ] = bheapInd
  h.lookup[ h.nodes[bheapInd].Id ] = aheapInd
  #Introducing temporary variables for the swapping.
  tempDist2::Float64 = h.nodes[ aheapInd ].dist2
  tempId::Int64 = h.nodes[ aheapInd ].Id
  #Swapping the entries of the corresponding node structs
  h.nodes[ aheapInd ].dist2 = h.nodes[ bheapInd ].dist2
  h.nodes[ aheapInd ].Id = h.nodes[ bheapInd ].Id 
  h.nodes[ bheapInd ].dist2 = tempDist2
  h.nodes[ bheapInd ].Id = tempId
end

#Compares the element with heap Index hInd to its two children and if it is smaller 
#than one of them, swaps it with the bigger one.
#The output is the the new hindex of the node, if a swap occured, and -1, otherwise.
function _moveDown!( h::Heap, hInd::Int64 )
  #Saving the squared distance of the node that is being swapped down.
  dist2::Float64 = h.nodes[ hInd ].dist2
  #If the both children exist:
  if 2 * hInd + 1 <= h.N 
    #If the left child is larger:
    if h.nodes[ 2 * hInd ].dist2 > h.nodes[ 2 * hInd + 1 ].dist2 
      #Check if the Child is larger than the parent:
      if h.nodes[2*hInd ].dist2 > dist2
        _swap!( h, hInd, 2 * hInd )
        return 2 * hInd
      else
        #No swap occuring:
        return -1
      end
    #If the left child is larger:
    else 
      #Check if the Child is larger than the parent:
      if h.nodes[ 2*hInd + 1 ].dist2 > dist2
        _swap!( h, hInd, 2 * hInd + 1 )
        return 2 * hInd + 1
      else
        #No swap occuring:
        return -1
      end
    end
  #If only one child exists: 
  elseif 2 * hInd <= h.N
    #Check if the Child is larger than the parent:
    if h.nodes[ 2 * hInd ].dist2 > dist2
      _swap!( h, hInd, 2 * hInd )
      return 2 * hInd
    end
  end
  #No swap occuring:
  return -1
end

function topId( h::Heap )
  return h.nodes[1].Id
end

function topDist2( h::Heap )
  return h.nodes[1].dist2
end

#Updates (decreases) an element of the heap and restores the heap property.
function update!( h::Heap, Id::Int64, newDist2::Float64 )
  tempInd::Int64 = h.lookup[ Id ]
  h.nodes[ tempInd ].dist2 = newDist2
  while ( tempInd <= h.N ) && ( tempInd >= 0 ) 
    tempInd = _moveDown!( h, tempInd )
  end
end

#############################################################################
#Introducing the "daycare", which keeps track of the descendants of every node
#############################################################################

mutable struct Daycare
  NParents::Int64
  NChildren::Int64
  NBuffer::Int64

  #This array gives for any parent node the corresponding ID
  IdParents::Array{Int64, 1}
  firstChildren::Array{Int64, 1}
  children::Array{Int64, 1}
  #This array has $N$ entries, that given node Id, the corresponding entry in the list 
  #of parents
  lookup::Array{Int64, 1} 
end

function newparent( dc::Daycare, IdParent::Int64 )
  dc.NParents = dc.NParents + 1
  dc.IdParents[ dc.NParents ] = IdParent
  dc.firstChildren[ dc.NParents ] = dc.NChildren + 1
  dc.firstChildren[ dc.NParents + 1 ] = dc.NChildren + 1
  dc.lookup[ IdParent ] = dc.NParents
end

function newson( dc::Daycare, IdSon::Int64 )
  if dc.NChildren >= dc.NBuffer
    dc.NBuffer = 2 * dc.NBuffer
    resize!( dc.children, dc.NBuffer ) 
  end
  dc.NChildren = dc.NChildren + 1
  dc.firstChildren[ dc.NParents + 1 ] = dc.firstChildren[ dc.NParents + 1 ] + 1
  dc.children[ dc.NChildren ] = IdSon
end

#############################################################################
#Implementation of the ordering algorithm
#############################################################################

#############################################################################
#Version without levels
#############################################################################
function _determineChildren!( h::Heap, 
                             dc::Daycare, 
                             parents::Array{Int64, 1}, 
                             x::Array{Float64, 2}, 
                             ID::Int64 )
#Function to determine the children of a node in the ordering
#Inputs:
# h:      
#   Heap to keep track of the elements and their distance to the points added
#   the ordering
# dc:
#   The "daycare" keeping track of the nodes and their children
# parents:
#   An array with N elements that in its i-th position contains the Id of the point
#   with Id i.
# x:      
#   A d times N array containing the positions of the data points
# ID:   
#   The ID of the point, the children of which are being determined

  #caching the Id of the father node
  father::Int64 = parents[ ID ]
  #caching the point
  y::Array{Float64, 1} = x[ :, ID ]
  #adding the a new node to the dc
  newparent( dc, ID )
  #cashing the ambient dimension
  d::Int64 = size( x, 1 )  
  #caching the distance of the newest elemnt
  pivotDist2::Float64 = topDist2( h )
  #We need to add the children to the new parent
  for index = dc.firstChildren[ dc.lookup[ father ] ]  : ( dc.firstChildren[ dc.lookup[ father ] + 1 ] - 1 )  
    j::Int64 = dc.children[ index ]
    dist2::Float64 = dist2eval( x, j , y, d )

    if dist2 <= pivotDist2
      if dist2 < h.nodes[ h.lookup[ j ] ].dist2 
        update!( h, j, dist2 )
      end
      newson( dc, j )
      if sqrt.( h.nodes[ h.lookup[ j ] ].dist2 ) .+ sqrt.( dist2 ) < sqrt.( pivotDist2 )
        parents[ j ] = ID
      end
    end
  end
end

function sortPoints( x::Array{Float64, 2}, boundDist::Array{Float64, 1} )
#Function to sort the points in a multicscale way
#Input:
# x:
#   A d times N array containing the measurement points
  @assert unique(x, 2) == x "Remove duplicate points first"
  Nbuffer0::Int64 = 10
  N::Int64 = size( x, 2 )
  parents::Array{Int64, 1} = Array{Int64, 1}(N)
  #Create heap:
  h::Heap = Heap( N, Array{HeapNode, 1}( N ), 1 : 1 : N )
  for i = 1 : N 
    h.nodes[ i ] = HeapNode( i, 1000000000000 )
  end
  #Create daycare
  dc::Daycare = Daycare( 0, 0, Nbuffer0, Array{Int64, 1}(N), Array{Int64, 1}(N+1),  Array{Int64, 1}(Nbuffer0), Array{Int64, 1}(N)) 
  #TODO remove ( only debugging )
  dc.IdParents .= 0
  dc.firstChildren .= 0
  dc.children .= 0
  dc.lookup .= 0
  parents .= 0
  #distances will contain as an i-th entry the distance upon elimination of the i-th
  #degree of freedom
  distances::Array{Float64} = zeros( N ) 


  #Setting the first element:
  #rootID::Int64 = randperm( N )[ 1 ]
  rootID::Int64 = 1;
  rootPoint::Array{Float64, 1} = x[ :, rootID ]
  parents .= rootID
  newparent( dc, rootID )
  for i = 1 : N
    newson( dc, i )
    distances[1] = max( distances[1], sqrt.( dist2eval( x, i, rootPoint, size( x, 1 ) ) ) )
    update!( h, i, dist2eval( x, i, rootPoint, size( x, 1 ) ) )
  end

  #updating the elements
  for iter = 2 : N
    pivotId::Int64 = topId( h )
    distances[ iter ] = sqrt( topDist2( h ) )
    _determineChildren!( h, dc, parents, x, pivotId )
  end

  #Output "P, revP, levels"
  #Adding a vector with the distances, at which the 

  return dc.IdParents, dc.lookup, distances
end

#Evaluates the distance of the point y in d-dimensional space to the jth point
#of the point cloud x
function dist2eval( x::Array{Float64, 2}, j::Int64, y::Array{Float64, 1}, d::Int64 )
  retVal::Float64 = 0.0;
  for i::Int64 = 1 : d
    retVal += ( x[ i, j ] - y[ i ] )^2
  end
  return retVal
end

#############################################################################
#Implementation of the sparsity pattern
#############################################################################

function _determineChildren!( revP::Array{Int64, 1}, distances::Array{Float64, 1}, dc::Daycare, parents::Array{Int64, 1}, x::Array{Float64, 2}, ID::Int64, rho::Float64 )
  #caches the father node
  father::Int64 = parents[ ID ]
  #caches the pivot point
  y::Array{Float64, 1} = x[ :, ID ]
  #Adds the new parent (row) to be added to the sparsity pattern
  newparent( dc, ID )
  #caches the ambient dimension
  d::Int64 = size( x, 1 )  
  #caches the distance coresponding to the pivot
  pivotDist::Float64 = distances[ revP[ ID ] ]
  for index = dc.firstChildren[ dc.lookup[ father ] ]  : ( dc.firstChildren[ dc.lookup[ father ] + 1 ] - 1 )  
    j::Int64 = dc.children[ index ]
    dist2::Float64 = dist2eval( x, j , y, d )

    if ( revP[ ID ] <= revP[ j ] ) && ( dist2 < rho^2 * ( pivotDist.^2 ) )
      newson( dc, j )
      if sqrt.( dist2 ) + rho * distances[ revP[ j ] ] <= rho * pivotDist  
        parents[ j ] = ID
      end
    end
  end
end



function sortPoints( x::Array{Float64, 2} )
#Function to sort the points in a multicscale way
#Input:
# x:
#   A d times N array containing the measurement points
  @assert unique(x, 2) == x "Remove duplicate points first"
  Nbuffer0::Int64 = 10
  N::Int64 = size( x, 2 )
  parents::Array{Int64, 1} = Array{Int64, 1}(N)
  #Create heap:
  h::Heap = Heap( N, Array{HeapNode, 1}( N ), 1 : 1 : N )
  for i = 1 : N 
    h.nodes[ i ] = HeapNode( i, 1000000000000 )
  end
  #Create daycare
  dc::Daycare = Daycare( 0, 0, Nbuffer0, Array{Int64, 1}(N), Array{Int64, 1}(N+1),  Array{Int64, 1}(Nbuffer0), Array{Int64, 1}(N)) 
  #TODO remove ( only debugging )
  dc.IdParents .= 0
  dc.firstChildren .= 0
  dc.children .= 0
  dc.lookup .= 0
  parents .= 0
  #distances will contain as an i-th entry the distance upon elimnation of the i-th
  #degree of freedom
  distances::Array{Float64} = zeros( N ) 


  #Setting the first element:
  #rootID::Int64 = randperm( N )[ 1 ]
  rootID::Int64 = 1;
  rootPoint::Array{Float64, 1} = x[ :, rootID ]
  parents .= rootID
  newparent( dc, rootID )
  for i = 1 : N
    newson( dc, i )
    distances[1] = max( distances[1], sqrt.( dist2eval( x, i, rootPoint, size( x, 1 ) ) ) )
    update!( h, i, dist2eval( x, i, rootPoint, size( x, 1 ) ) )
  end

  #updating the elements
  for iter = 2 : N
    pivotId::Int64 = topId( h )
    distances[ iter ] = sqrt( topDist2( h ) )
    _determineChildren!( h, dc, parents, x, pivotId )
  end

  #Output "P, revP, levels"
  #Adding a vector with the distances, at which the 

  return dc.IdParents, dc.lookup, distances
end

#############################################################################
#Implementation of the sparsity pattern
#############################################################################

function sparsityPattern( x::Array{Float64, 2}, P::Array{Int64, 1}, revP::Array{Int64, 1}, distances::Array{Float64, 1}, rho::Float64 )
  d::Int64 = size( x, 1 )
  N::Int64 = size( x, 2 )
  Nbuffer0::Int64 = 10

  parents::Array{Int64, 1} = Array{Int64, 1}(N)
  dc::Daycare = Daycare( 0, 0, Nbuffer0, Array{Int64, 1}(N), Array{Int64, 1}(N+1),  Array{Int64, 1}(Nbuffer0), Array{Int64, 1}(N)) 

  dc.IdParents .= 0
  dc.firstChildren .= 0
  dc.children .= 0
  dc.lookup .= 0
  parents .= 0

  parents .= P[ 1 ]
  newparent( dc, P[ 1 ] )
  for i = 1 : N
    newson( dc, i )
  end

  for k = 2 : N
    _determineChildren!( revP, distances, dc, parents, x, P[ k ], rho)
  end
  
  #Sort the jind according to their ordering:
  for i = 1 : N
    dc.children[ dc.firstChildren[ i ] : ( dc.firstChildren[ i + 1 ] - 1 ) ] = P[ sort!( revP[ dc.children[ dc.firstChildren[ i ] : ( dc.firstChildren[ i + 1 ] - 1 ) ] ] ) ]
  end

  return dc.firstChildren, revP[ dc.children[ 1 : dc.NChildren ] ]
end

