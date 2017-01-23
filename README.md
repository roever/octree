# Sparse 3D-Array
Sparse 3d array implemented as octree

The octree is implemented inside the class itself and not via lots of cross linked object. All nodes of the tree are stored within 2 vectors
(one for the inside nodes, one for the leafs) increasing cache locality and hopefully speed (non testet though)

The array will automatically grow up to the size required for the coordinates of the stored elements.

Types for the dimensions (e.g., how big is the adressable volume), types for storage (e.g. how many none-null nodes are maximally inside the tree)
can be specified via template parameters



