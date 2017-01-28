#pragma once

#include <cstdint>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <limits>
#include <stdexcept>
#include <array>

#if __has_include(<boost/qvm/vec_access.hpp>)
#include <boost/qvm/vec_access.hpp>
#endif

template <class I> struct fast_type_for_int;

template <> struct fast_type_for_int<uint8_t>  { using type=uint_fast8_t;  };
template <> struct fast_type_for_int<int8_t>   { using type=int_fast8_t;   };
template <> struct fast_type_for_int<uint16_t> { using type=uint_fast16_t; };
template <> struct fast_type_for_int<int16_t>  { using type=int_fast16_t;  };
template <> struct fast_type_for_int<uint32_t> { using type=uint_fast32_t; };
template <> struct fast_type_for_int<int32_t>  { using type=int_fast32_t;  };
template <> struct fast_type_for_int<uint64_t> { using type=uint_fast64_t; };
template <> struct fast_type_for_int<int64_t>  { using type=int_fast64_t;  };

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
/// \tparam D: type to use for accessing the dimensions this limits the maximal dimensions
///            you can access. It must be a signed type
/// \tparam I: index type used as pointers for leafs and nodes of the tree. This type and B
///            defines how many non-default values can be stored. D should be an unsigned type
///            for maximal range. Assuming M is the maximal value storable in I you can
///            store between M/8 and M nun default values in the structure
///            if things run over, an std::runtime_error is thrown
/// \tparam II: type used internally to represent I, this type should be at least as big as I
///             but might be bigger, if that type is faster to handle than I (e.g. use the fast_uintxx_t
///             types) The deault for this template parameter automatically gets the fast_int_t type
///             suitable for I.
template<class T, int B=2, class D=int_fast32_t, class I=uint32_t, class II=typename fast_type_for_int<I>::type>
class Sparse3DArray
{
  private:
    // the datatype used to indes into the nodes and leafs vectors
    using index_type=I;

    // the inside nodes of the octree always 8 consecutive indices are the children of one node
    // the first node with index 0 is always the root node, which will never be freed
    // making the minimal size of the tree B (from -B <= d < B is accessible)
    std::vector<I> nodes;
    // the leaf nodes of the octree always B*B*B consecutive entries create one leaf node,
    // entry 0 is the empty value
    std::vector<T> leafs;
    // size of the tree, area that is accessible in the 3 dimensions (from -size <= d < size)
    D size;
    // index of the first empty leaf or zero empy nodes are handled by having a sinly linked
    // list of nodes. for nodes the first index in there points to the next
    // empty node, for leafs a bit of casting is required. That is the reason why
    // we need B to be big enough
    II firstEmptyLeaf;
    // index of the first empty node or zero
    II firstEmptyNode;

    // function to allocate a leaf or a node, the returned value will always be initialiazed
    // to init
    template <class C>
    II alloc(C & c, II & firstEmpty, II n, typename C::value_type init)
    {
      if (firstEmpty)
      {
        // we have a node in the empty list, so take it from there
        II result = firstEmpty;
        firstEmpty = *reinterpret_cast<I*>(&c[result]);

        // reset it to init
        for (II i = 0; i < n; i++)
          c[result+i] = init;

        return result;
      }
      else
      {
        // no empty node available, resize the vector, but first check, if we can still
        // access the new elements
        size_t result = c.size();
        if (result >= std::numeric_limits<I>::max()-n) throw std::runtime_error("Too many nodes in octree");
        c.resize(result+n, init);
        return result;
      }
    }

    // free a node by adding it to the empty list
    template <class C>
    void free(II i, C & c, II & firstEmpty) noexcept
    {
      *reinterpret_cast<I*>(&c[i]) = firstEmpty;
      firstEmpty = i;
    }

    // functions to allocate and free the nodes and leafs
    II allocLeaf() { return alloc(leafs, firstEmptyLeaf, B*B*B, T()); }
    II allocNode() { return alloc(nodes, firstEmptyNode, 8, 0); }

    void freeLeaf(II i) noexcept { free(i, leafs, firstEmptyLeaf); }
    void freeNode(II i) noexcept { free(i, nodes, firstEmptyNode); }

    // check, if (x, y, z) is inside the current size of the tree
    bool inside(D x, D y, D z) const noexcept
    {
      return -size <= x && x < size && -size <= y && y < size && -size <= z && z < size;
    }

    // calculate the index into the leaf node, assuming the given
    // coordinates are centered around zero
    constexpr II leafIdx(D x, D y, D z) const noexcept
    {
      return ((z+B/2)*B+(y+B/2))*B+(x+B/2);
    }

    // calculate index into a node of the octree
    // size is assumed to be the size of the current level of the octree
    // and is automatically halved for the next level
    // x, y and z are also updated to the coordinates that need to be used
    // when recursively calling the next function on the subtree
    // the returned index is the one where the subtree can be found
    constexpr II nodeIdx(II idx, D & x, D & y, D & z, D & size) const noexcept
    {
      size /= 2;
      if (z >= 0) { idx += 4; z -= size; } else { z += size; }
      if (y >= 0) { idx += 2; y -= size; } else { y += size; }
      if (x >= 0) { idx += 1; x -= size; } else { x += size; }

      return idx;
    }

    // recursive get-function, assumes that (x, y, z) is inside of the tree
    const T & recGet(D x, D y, D z, II idx, D size) const noexcept
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
    void recSet(D x, D y, D z, II idx, D size, const T & val)
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
          II i = (size >= B) ? allocNode() : allocLeaf();
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
    bool recReset(D x, D y, D z, II idx, D size) noexcept
    {
      if (size >= B)
      {
        II idx2 = nodeIdx(idx, x, y, z, size);

        if (nodes[idx2])
        {
          // reset the child
          if (recReset(x, y, z, nodes[idx2], size))
          {
            // child has been cleared... so remove the pointer and check
            // if the complete node is now empty
            nodes[idx2] = 0;
            for (II i = 0; i < 8; i++)
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
        II idx2 = idx + leafIdx(x, y, z);
        if (leafs[idx2] != T())
        {
          // the value is not clear, so clear is and check if the leaf
          // is now completely empty, if so free it and tell the caller
          leafs[idx2] = T();
          for (II i = 0; i < B*B*B; i++)
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
      for (II i = 0; i < 8; i++)
        if (nodes[i])
        {
          II n = allocNode();
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
        for (II i = 0; i < 8; i++)
          if (nodes[i])
            for (II j = 0; j < 8; j++)
              if ((i^7) != j)
                if (nodes[nodes[i]+j])
                  return;

        // outer layers are empty, drop one level
        for (II i = 0; i < 8; i++)
        {
          if (nodes[i])
          {
            II n = nodes[i];
            nodes[i] = nodes[nodes[i]+(i^7)];
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
    void recForEach(F f, D x, D y, D z, II idx, D size) const noexcept(std::is_nothrow_callable(F::operator()))
#else
    void recForEach(F f, D x, D y, D z, II idx, D size) const
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
        II i = 0;

        for (D iz = -B/2; iz < B/2; iz++)
          for (D iy = -B/2; iy < B/2; iy++)
            for (D ix = -B/2; ix < B/2; ix++)
            {
              if (leafs[idx+i] != T())
                f(x+ix, y+iy, z+iz, leafs[idx+i]);
              i++;
            }
      }
    }

    // recursively calculate the bounding box, it is done by calculating the box for all
    // the children and then assemble the big box out of that information
    // edges is a bitmask of the edges that need to be found,
    // cx, cy, zy is the center of the octree, size the size of the cube
    // the function will change the bounding box in sz according to what is found
    void recCalcBoundingBox(std::array<D, 6> & sz, D cx, D cy, D cz, II idx, D size) const
    {
      // allright, this warrants some more explanation, the recursive function enlarges the
      // given bounding box, if needed, so you need to instantiate it properly

      // before we start looking though, we have a look at the already found edges, we only
      // start our search, if we can actually improve them. If the currently found edge is
      // already beyond the size of the current tree, there is no need to do anything
      if (   (sz[0] <= cx-size) && (sz[1] >= cx+size-1)
          && (sz[2] <= cy-size) && (sz[3] >= cy+size-1)
          && (sz[4] <= cz-size) && (sz[5] >= cz+size-1)
         )
      {
        // if no edge is requested, just do nothing
        return;
      }
      else if (size >= B)
      {
        size /= 2;

        // calculate the bounding box of the children.
        // the order of the child sears is done in such a way that we can hope to get close
        // to the final box very quickly
        if (nodes[idx+0]) { recCalcBoundingBox(sz, cx-size, cy-size, cz-size, nodes[idx+0], size); }
        if (nodes[idx+7]) { recCalcBoundingBox(sz, cx+size, cy+size, cz+size, nodes[idx+7], size); }

        if (nodes[idx+1]) { recCalcBoundingBox(sz, cx+size, cy-size, cz-size, nodes[idx+1], size); }
        if (nodes[idx+6]) { recCalcBoundingBox(sz, cx-size, cy+size, cz+size, nodes[idx+6], size); }

        if (nodes[idx+2]) { recCalcBoundingBox(sz, cx-size, cy+size, cz-size, nodes[idx+2], size); }
        if (nodes[idx+5]) { recCalcBoundingBox(sz, cx+size, cy-size, cz+size, nodes[idx+5], size); }

        if (nodes[idx+3]) { recCalcBoundingBox(sz, cx+size, cy+size, cz-size, nodes[idx+3], size); }
        if (nodes[idx+4]) { recCalcBoundingBox(sz, cx-size, cy-size, cz+size, nodes[idx+4], size); }
      }
      else
      {
        // for the leaf iterate over the content
        II i = 0;

        for (D iz = -B/2; iz < B/2; iz++)
          for (D iy = -B/2; iy < B/2; iy++)
            for (D ix = -B/2; ix < B/2; ix++)
            {
              if (leafs[idx+i] != T())
              {
                sz[0] = std::min(sz[0], cx+ix);
                sz[1] = std::max(sz[1], cx+ix);
                sz[2] = std::min(sz[2], cy+iy);
                sz[3] = std::max(sz[3], cy+iy);
                sz[4] = std::min(sz[4], cz+iz);
                sz[5] = std::max(sz[5], cz+iz);
              }
              i++;
            }
      }
    }

    // free all data of the tree, resetting in the end everything
    // to default value
    void recClear(II idx, D size) noexcept
    {
      if (size > B)
      {
        size /= 2;
        for (II i = 0; i < 8; i++)
          if (nodes[idx+i])
          {
            recClear(nodes[idx+i], size);
            freeNode(nodes[idx+i]);
            nodes[idx+i] = 0;
          }
      }
      else
      {
        for (II i = 0; i < 8; i++)
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
      static_assert(sizeof(T)*B*B*B >= sizeof(I), "B must be big enough to fit an index into a leaf");
      nodes.resize(8);
      leafs.resize(1);
    }

    /// get the value of an entry at the given position
    /// \param x x-coordinate to get
    /// \param y y-coordinate to get
    /// \param z z-coordinate to get
    const T & get(D x, D y, D z) const noexcept
    {
      if (inside(x, y, z))
        return recGet(x, y, z, 0, size);
      else
        return leafs[0];
    }

    /// set a value, a copy is stored at the given position
    /// when the structure is "full", meaning it is not possible to
    /// add more adressable nodes it will throw a std::runtime_exception
    void set(D x, D y, D z, const T & val)
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

#if __has_include(<boost/qvm/vec_access.hpp>)
    /// like get but you can use any boost qvm valid vector to specify
    /// the position... TODO if non-integer members in vector round them
    template <class V>
    const T & get(const V & v) const noexcept
    {
      return get(boost::qvm::X(v), boost::qvm::Y(v), boost::qvm::Z(v));
    }

    /// like set but you can use any boost qvm valid vector to specify
    /// the position... TODO if non-integer members in vector round them
    template <class V>
    void set(const V & v, const T & val) const noexcept
    {
      return set(boost::qvm::X(v), boost::qvm::Y(v), boost::qvm::Z(v), val);
    }
#endif

    /// call f for each non default value within the space, f must
    /// have 4 arguments, (D x, D y, D z, const T & t) or
    /// (D x, D y, D z, T t), if your T is simple to copy.
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
    void swap(Sparse3DArray & a) noexcept(noexcept(std::swap<std::vector<T>>) && noexcept(std::swap<std::vector<I>>))
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
    bool calcBoundingBox(D & xmin, D & ymin, D & zmin, D & xmax, D & ymax, D & zmax) const noexcept
    {
      // initialize the bounding box so that it is empty, the recursive
      // function below will move the boundary according to what is found
      std::array<D, 6> sz = { { size, -size, size, -size, size, -size } };

      recCalcBoundingBox(sz, 0, 0, 0, 0, size);

      xmin = sz[0];
      xmax = sz[1];
      ymin = sz[2];
      ymax = sz[3];
      zmin = sz[4];
      zmax = sz[5];

      return xmin <= xmax;
    }

    /// clear everything back to default constructor values. The allocated memory will be kepts though
    void clear() noexcept
    {
      recClear(0, size);
      merge();
    }
};

