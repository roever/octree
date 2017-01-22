#pragma once

#include <cstdint>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <limits>
#include <stdexcept>

/// Implement a 3d-array. The array is self-growing and can potentially
/// contain all int-adressable 3d-coordinates. The array contains elements
/// of type T. The default constructed T() is special, all cells never
/// accessed have this value and don't use up any memory. All other values
/// do. 
/// \note Internally it is implemented as an Octree
///
/// \tparam T: datatype to put into this tree, keep it a simple datatype, because internally
///            there are some places where the value is copied
/// \tparam B: merged base level, size of one dimension of the base level, so one cube
///            at the base is always B*B*B elements of T, must be even
///            B must be big enough so that the resulting B*B*B instances of T can contain
///            one instance of I
/// \tparam I: index type used as pointers for leafs and nodes of the tree, if the type is too
///            small, the array can not contain many non-default values
template<class T, int B=2, class I=uint32_t>
class Sparse3DArray
{
  private:
    // the datatype used to indes into the nodes and leafs vectors
    using index_type=I;

    // the inside nodes of the octree always 8 consecutive indices are the children of one node
    // the first node with index 0 is always the root node, which will never be freed
    // making the minimal size of the tree B (from -B <= d < B is accessible)
    std::vector<index_type> nodes;  
    // the leaf nodes of the octree always B*B*B consecutive entries create one leaf node, 
    // entry 0 is the empty value
    std::vector<T> leafs;  
    // size of the tree, area that is accessible in the 3 dimensions (from -size <= d < size)
    int size;       
    // index of the first empty leaf or zero
    int firstEmptyLeaf;
    // index of the first empty node or zero
    int firstEmptyNode;

    // function to allocate a leaf or a node, the returned value will always be initialiazed
    // to init
    template <class C>
    int alloc(C & c, int & firstEmpty, int n, typename C::value_type init)
    {
      if (firstEmpty)
      {
        // we have a node in the empty list, so take it from there
        int result = firstEmpty;
        firstEmpty = *(index_type*)&c[result];

        // reset it to init
        for (int i = 0; i < n; i++)
          c[result+i] = init;

        return result;
      }
      else
      {
        // no empty node available, resize the vector, but first check, if we can still
        // access the new elements
        size_t result = c.size();
        if (result >= std::numeric_limits<index_type>::max()-n) throw std::runtime_error("Too many nodes in octree");
        c.resize(result+n, init);
        return result;
      }
    }

    // free a node by adding it to the empty list
    template <class C>
    void free(int i, C & c, int & firstEmpty) noexcept
    {
      *(index_type*)&c[i] = firstEmpty;
      firstEmpty = i;
    }

    // functions to allocate and free the nodes and leafs
    int allocLeaf() { return alloc(leafs, firstEmptyLeaf, B*B*B, T()); }
    int allocNode() { return alloc(nodes, firstEmptyNode, 8, 0); }
    
    void freeLeaf(int i) noexcept { free(i, leafs, firstEmptyLeaf); }
    void freeNode(int i) noexcept { free(i, nodes, firstEmptyNode); }

    // check, if (x, y, z) is inside the current size of the tree
    bool inside(int x, int y, int z) const noexcept
    {
      return -size <= x && x < size && -size <= y && y < size && -size <= z && z < size;
    }

    // calculate the index into the leaf node, assuming the given
    // coordinates are centered around zero
    constexpr int leafIdx(int x, int y, int z) const noexcept
    {
      return ((z+B/2)*B+(y+B/2))*B+(x+B/2);
    }

    // calculate index into a node of the octree
    // size is assumed to be the size of the current level of the octree
    // and is automatically halved for the next level
    // x, y and z are also updated to the coordinates that need to be used
    // when recursively calling the next function on the subtree
    // the returned index is the one where the subtree can be found
    constexpr int nodeIdx(int idx, int & x, int & y, int & z, int & size) const noexcept
    {
      size /= 2;
      if (z >= 0) { idx += 4; z -= size; } else { z += size; }
      if (y >= 0) { idx += 2; y -= size; } else { y += size; }
      if (x >= 0) { idx += 1; x -= size; } else { x += size; }

      return idx;
    }

    // recursive get-function, assumes that (x, y, z) is inside of the tree
    const T & recGet(int x, int y, int z, int idx, int size) const noexcept
    {
      if (size >= B)
      {
        idx = nodeIdx(idx, x, y, z, size);

        if (nodes[idx])
          return recGet(x, y, z, nodes[idx], size);
        else
          return leafs[0];
      }
      else
      {
        return leafs[idx+leafIdx(x, y, z)];
      }
    }

    // sets a value, creating the required nodes on the way, assumes
    // that (x, y, z) is inside the tree and val is not the default
    // constructible value
    void recSet(int x, int y, int z, int idx, int size, const T & val)
    {
      if (size >= B)
      {
        idx = nodeIdx(idx, x, y, z, size);

        if (!nodes[idx])
        {
          // big oups here, if we assign directly without going over i, we get a race
          // because allocNode might resize the nodes container and the left side of the
          // assignment might be invalidated leading to an access into heap memory that had
          // been freed...
          int i = (size >= B) ? allocNode() : allocLeaf();
          nodes[idx] = i;
        }

        recSet(x, y, z, nodes[idx], size, val);
      }
      else
      {
        leafs[idx+leafIdx(x, y, z)] = val;
      }
    }

    // resets a field back to the default value. Assumes (x, y, z) is inside the
    // tree. Frees empty nodes on the way up from the recursion.
    // As we don't want to free the root node, each recursion check only the children
    // and frees them, but not itself, the function returns true, 
    // if all elements of the child are empty and thus the child can
    // be freed
    bool recReset(int x, int y, int z, int idx, int size) noexcept
    {
      if (size >= B)
      {
        int idx2 = nodeIdx(idx, x, y, z, size);

        if (nodes[idx2])
        {
          // reset the child
          if (recReset(x, y, z, nodes[idx2], size))
          {
            // child has been cleared... so remove the pointer and check
            // if the complete node is now empty
            nodes[idx2] = 0;
            for (int i = 0; i < 8; i++)
              if (nodes[idx+i])
                return false;

            // if so, and the node is not root, free it
            if (idx)
              freeNode(idx);

            // and tell the called that the node was free and has been
            // cleared
            return true;
          }

          // either nothing has been changed in child, or it is still not
          // empty, anyways, we don't need to do anything special
          return false;
        }
        else
        {
          // nothing has been changed, node is already clear, 
          // so no need to check anyting
          return false;
        }
      }
      else
      {
        int idx2 = idx + leafIdx(x, y, z);
        if (leafs[idx2] != T())
        {
          // the value is not clear, so clear is and check if the leaf
          // is now completely empty, if so free it and tell the caller
          leafs[idx2] = T();
          for (int i = 0; i < B*B*B; i++)
            if (leafs[idx+i] != T())
              return false;

          freeLeaf(idx);
          return true;
        }

        // leaf is already clear
        return false;
      }
    }

    // increase level of the ocreee by one, doubling the accessible 
    // size. The content will stay where is is, it is just one
    // level deeper in the tree
    void split(void)
    {
      for (int i = 0; i < 8; i++)
        if (nodes[i])
        {
          int n = allocNode();
          // the xor 7 will invert the lowest 3 bits in practice giving
          // us the index of the opposite subtree of the node
          nodes[n+(i^7)] = nodes[i];
          nodes[i] = n;
        }

      size *= 2;
    }

    // try to remove as many levels as possible from the octree without
    // loosing information
    void merge(void) noexcept
    {
      while (size > B)
      {
        // first check, whether the outer layer subtrees are all empty
        // if not leave the function, there is nothing we can do
        for (int i = 0; i < 8; i++)
          if (nodes[i])
            for (int j = 0; j < 8; j++)
              if ((i^7) != j)
                if (nodes[nodes[i]+j])
                  return;

        // outer layers are empty, drop one level
        for (int i = 0; i < 8; i++)
        {
          if (nodes[i])
          {
            int n = nodes[i];
            nodes[i] = nodes[nodes[i]+i^7];
            freeNode(n);
          }
        }

        // accessible size has been halfed
        size /= 2;
      }
    }

    // iterate over the tree and call f for all non default entries
    // x, y, z point to the current center of the subtree
    template<class F>
#ifdef USE_CPP_17
    void recForEach(F f, int x, int y, int z, int idx, int size) const noexcept(std::is_nothrow_callable(F::operator()))
#else
    void recForEach(F f, int x, int y, int z, int idx, int size) const 
#endif
    {
      if (size >= B)
      {
        size /= 2;
        // check all 8 children and recursively call function, if they are occupied
        if (nodes[idx+0]) recForEach(f, x-size, y-size, z-size, nodes[idx+0], size);
        if (nodes[idx+1]) recForEach(f, x+size, y-size, z-size, nodes[idx+1], size);
        if (nodes[idx+2]) recForEach(f, x-size, y+size, z-size, nodes[idx+2], size);
        if (nodes[idx+3]) recForEach(f, x+size, y+size, z-size, nodes[idx+3], size);
        if (nodes[idx+4]) recForEach(f, x-size, y-size, z+size, nodes[idx+4], size);
        if (nodes[idx+5]) recForEach(f, x+size, y-size, z+size, nodes[idx+5], size);
        if (nodes[idx+6]) recForEach(f, x-size, y+size, z+size, nodes[idx+6], size);
        if (nodes[idx+7]) recForEach(f, x+size, y+size, z+size, nodes[idx+7], size);
      }
      else
      {
        // leaf level reached, iterate over the content of the leaf
        for (int i = 0; i < B*B*B; i++)
          if (leafs[idx+i] != T())
            f(x+i%B-B/2, y+(i/B)%B-B/2, z+(i/(B*B))-B/2, leafs[idx+i]);
      }
    }

    // recursively calculate the bounding box, it is done by calculating the box for all
    // the children and then assemble the big box out of that information
    // return true, when bounding box is valid
    int recCalcBoundingBox(int & xmin, int & ymin, int & zmin, int & xmax, int & ymax, int & zmax, int idx, int size) const
    {
      bool result = false;

      if (size >= B)
      {
        size /= 2;

        // go over all children
        for (int i = 0; i < 8; i++)
          if (nodes[idx+i])
          {
            int x1, y1, z1, x2, y2, z2;

            // calculate their bounding boxes
            if (recCalcBoundingBox(x1, y1, z1, x2, y2, z2, nodes[idx+i], size))
            {
              // move them to the right position
              if ((i & 1) == 0) { x1 -= size; x2 -= size; } else { x1 += size; x2 += size; }
              if ((i & 2) == 0) { y1 -= size; y2 -= size; } else { y1 += size; y2 += size; }
              if ((i & 4) == 0) { z1 -= size; z2 -= size; } else { z1 += size; z2 += size; }

              // either take over the values, if we don't have our own box, yet
              // of make our box larger to include child
              if (!result)
              {
                xmin = x1; ymin = y1; zmin = z1;
                xmax = x2; ymax = y2; zmax = z2;
                result = true;
              }
              else
              {
                xmin = std::min(xmin, x1);
                ymin = std::min(ymin, y1);
                zmin = std::min(zmin, z1);

                xmax = std::max(xmax, x2);
                ymax = std::max(ymax, y2);
                zmax = std::max(zmax, z2);
              }
            }
          }
      }
      else
      {
        // for the leaf iterate over the content
        xmin = B; xmax = -B;
        ymin = B; ymax = -B;
        zmin = B; zmax = -B;

        int i = 0;

        for (int iz = -B/2; iz < B/2; iz++)
          for (int iy = -B/2; iy < B/2; iy++)
            for (int ix = -B/2; ix < B/2; ix++)
            {
              if (leafs[idx+i] != T())
              {
                xmin = std::min(xmin, ix);
                xmax = std::max(xmax, ix);
                ymin = std::min(ymin, iy);
                ymax = std::max(ymax, iy);
                zmin = std::min(zmin, iz);
                zmax = std::max(zmax, iz);

                result = true;
              }
              i++;
            }
      }

      return result;
    }
      
    // free all data of the tree, resetting in the end everything
    // to default value
    void recClear(int idx, int size) noexcept
    {
      if (size > B)
      {
        size /= 2;
        for (int i = 0; i < 8; i++)
          if (nodes[idx+i])
          {
            recClear(nodes[idx+i], size);
            freeNode(nodes[idx+i]);
            nodes[idx+i] = 0;
          }
      }
      else
      {
        for (int i = 0; i < 8; i++)
          if (nodes[idx+i])
          {
            freeLeaf(nodes[idx+i]);
            nodes[idx+i] = 0;
          }
      }
    }


  public:

    /// create an empty array, all values in the 3d-space are default initialized
    Sparse3DArray() : size(B), firstEmptyLeaf(0), firstEmptyNode(0)
    {
      static_assert((B % 2) == 0, "B must be even");
      static_assert(sizeof(T)*B*B*B >= sizeof(index_type), "B must be big enough to fit an index into a leaf");
      nodes.resize(8);
      leafs.resize(1);
    }

    /// get the value of an entry at the given position
    /// \param x x-coordinate to get
    /// \param y y-coordinate to get
    /// \param z z-coordinate to get
    const T & get(int x, int y, int z) const noexcept
    {
      if (inside(x, y, z))
        return recGet(x, y, z, 0, size);
      else
        return leafs[0];
    }

    /// set a value, a copy is stored at the given position
    /// when the structure is "full", meaning it is not possible to
    /// add more adressable nodes it will throw a std::runtime_exception
    void set(int x, int y, int z, const T & val) 
    {
      if (val != T())
      {
        while (!inside(x, y, z)) { split(); }
        recSet(x, y, z, 0, size, val);
      }
      else if (inside(x, y, z))
      {
        // the user wants to reset a value to default, so
        // call the reset function and when it returns
        // with true, meaning there is something that was cleared
        // in the root node we try to merge. This is conservative
        // because we might be able to merge sooner, but that is not
        // that easy to find out
        if (recReset(x, y, z, 0, size))
          merge();
      }
      else
      {
        // outside of tree, no need to clear things, already clear
      }
    }

    /// call f for each non default value within the space, f must
    /// have 4 arguments, (int x, int y, int z, const T & t) or
    /// (int x, int y, int z, T t), if your T is simple to copy.
    /// You can not make any modifications while running for_each
    /// if you do modifications while this is running you might end up
    /// called for some old and some new values. You better collect
    /// all the changes you want to make and apply them afterwards
    /// \tparam F callable type
    /// \param f callable object
    template<class F>
#ifdef USE_CPP_17
    void for_each(F f) const noexcept(std::is_nothrow_callable(F::operator()))
#else
    void for_each(F f) const 
#endif
    {
      recForEach(f, 0, 0, 0, 0, size);
    }

    /// swap two array
    /// \param a the array to swap with
    void swap(Sparse3DArray & a) noexcept(noexcept(std::swap<std::vector<T>>) && noexcept(std::swap<std::vector<index_type>>))
    {
      std::swap(nodes, a.nodes);
      std::swap(leafs, a.leafs);
      std::swap(size, a.size);
      std::swap(firstEmptyLeaf, a.firstEmptyLeaf);
      std::swap(firstEmptyNode, a.firstEmptyNode);
    }

    /// calculate the bounding box
    /// \param xmin minimal x-coordinate found
    /// \param ymin minimal y-coordinate found
    /// \param zmin minimal z-coordinate found
    /// \param xmax maximal x-coordinate found
    /// \param ymax maximal y-coordinate found
    /// \param zmax maximal z-coordinate found
    /// \return true, when the array contains data and the bounding box is value, false if the array is empty
    bool calcBoundingBox(int & xmin, int & ymin, int & zmin, int & xmax, int & ymax, int & zmax) const noexcept
    {
      return recCalcBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax, 0, size);
    }

    /// clear everything back to default constructor values. The allocated memory will be kepts though
    void clear() noexcept
    {
      recClear(0, size);
      merge();
    }
};

