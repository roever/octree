#include "octree.hpp"

#define NN 500

#include <cstdint>
#include <cstdio>
#include <ctime>

int main()
{
  Sparse3DArray<std::int16_t, 2> array2;
  Sparse3DArray<std::int16_t, 4> array4;
  Sparse3DArray<std::int16_t, 6> array6;
  Sparse3DArray<std::int16_t, 8> array8;

  array2.set(4,5,6,1);
  array2.set(102,321,678,4);
  array4.set(4,5,6,1);
  array4.set(102,321,678,4);
  array6.set(4,5,6,1);
  array6.set(102,321,678,4);
  array8.set(4,5,6,1);
  array8.set(102,321,678,4);

  array2.for_each([](int x, int y, int z, int val) { printf("%i %i %i = %i\n", x, y, z, val); });
  array4.for_each([](int x, int y, int z, int val) { printf("%i %i %i = %i\n", x, y, z, val); });
  array6.for_each([](int x, int y, int z, int val) { printf("%i %i %i = %i\n", x, y, z, val); });
  array8.for_each([](int x, int y, int z, int val) { printf("%i %i %i = %i\n", x, y, z, val); });

  const auto a2 = array2;

  a2.for_each([](int x, int y, int z, int val) { printf("%i %i %i = %i\n", x, y, z, val); });
  array2.for_each([](int x, int y, int z, const int & val) { printf("%i %i %i = %i\n", x, y, z, val); });

  Sparse3DArray<std::int16_t, 2> array2_2;

  array2_2.swap(array2);
  std::swap(array2, array2_2);


  int x1, x2, y1, y2, z1, z2;

  array2.calcBoundingBox(x1, y1, z1, x2, y2, z2); printf("bb %i-%i, %i-%i, %i-%i\n", x1, x2, y1, y2, z1, z2);
  array4.calcBoundingBox(x1, y1, z1, x2, y2, z2); printf("bb %i-%i, %i-%i, %i-%i\n", x1, x2, y1, y2, z1, z2);
  array6.calcBoundingBox(x1, y1, z1, x2, y2, z2); printf("bb %i-%i, %i-%i, %i-%i\n", x1, x2, y1, y2, z1, z2);
  array8.calcBoundingBox(x1, y1, z1, x2, y2, z2); printf("bb %i-%i, %i-%i, %i-%i\n", x1, x2, y1, y2, z1, z2);

  int loops = 1;

  std::vector<int> c;

#if 0
  while (true)
  {
    if (loops % 100 == 0)
      printf("loop: %i\n", loops);
    loops++;

    bool bad = false;

    do {

      c.clear();
      for (int i = 0; i < 1000; i++)
        c.push_back(rand()%2000-1000);

      bad = false;
      for (int i = 1; i < NN-4; i++)
        for (int j = 0; j < i; j++)
        {
          if (c[i] == c[j] && c[i+1] == c[j+1] && c[i+2] == c[j+2] && c[i+3] != c[j+3])
          {
            bad = true;
            i = NN;
            break;
          }
        }

    } while (bad);

    for (int i = 0; i < NN-4; i++)
    {
      //printf("s %i %i %i  %i\n", c[i], c[i+1], c[i+2], c[i+3]);
      array2.set(c[i], c[i+1], c[i+2], c[i+3]);
      array4.set(c[i], c[i+1], c[i+2], c[i+3]);
      array6.set(c[i], c[i+1], c[i+2], c[i+3]);
      array8.set(c[i], c[i+1], c[i+2], c[i+3]);
    }

    for (int i = 0; i < NN-4; i++)
    {
      //printf("t %i %i %i  %i\n", c[i], c[i+1], c[i+2], c[i+3]);

      if (array2.get(c[i], c[i+1], c[i+2]) != c[i+3]) exit(1);
      if (array4.get(c[i], c[i+1], c[i+2]) != c[i+3]) exit(1);
      if (array6.get(c[i], c[i+1], c[i+2]) != c[i+3]) exit(1);
      if (array8.get(c[i], c[i+1], c[i+2]) != c[i+3]) exit(1);
    }

    for (int i = 0; i < NN-4; i++)
    {
      //printf("s %i %i %i  %i\n", c[i], c[i+1], c[i+2], c[i+3]);
      array2.set(c[i], c[i+1], c[i+2], 0);
      array4.set(c[i], c[i+1], c[i+2], 0);
      array6.set(c[i], c[i+1], c[i+2], 0);
      array8.set(c[i], c[i+1], c[i+2], 0);
    }

    for (int i = 0; i < NN-4; i++)
    {
      if (array2.get(c[i], c[i+1], c[i+2]) != 0) exit(1);
      if (array4.get(c[i], c[i+1], c[i+2]) != 0) exit(1);
      if (array6.get(c[i], c[i+1], c[i+2]) != 0) exit(1);
      if (array8.get(c[i], c[i+1], c[i+2]) != 0) exit(1);
    }
  }
#endif

  time_t start = time(NULL);

  while (true)
  {
    if (loops % 1000 == 0)
      printf("loop: %i %f\n", loops, 1000.0*(time(NULL)-start)/loops);
    loops++;

    bool bad = false;

    do {

      c.clear();
      for (int i = 0; i < NN; i++)
        c.push_back(rand()%200-100);

      bad = false;
      for (int i = 1; i < NN-4; i++)
        for (int j = 0; j < i; j++)
        {
          if (c[i] == c[j] && c[i+1] == c[j+1] && c[i+2] == c[j+2] && c[i+3] != c[j+3])
          {
            bad = true;
            i = NN;
            break;
          }
        }

    } while (bad);

    for (int i = 0; i < NN-4; i++)
    {
      array2.set(c[i], c[i+1], c[i+2], c[i+3]);
    }

    for (int i = 0; i < NN-4; i++)
    {
      if (array2.get(c[i], c[i+1], c[i+2]) != c[i+3]) exit(1);
    }

    for (int i = 0; i < NN-4; i++)
    {
      array2.set(c[i], c[i+1], c[i+2], 0);
    }

    for (int i = 0; i < NN-4; i++)
    {
      if (array2.get(c[i], c[i+1], c[i+2]) != 0) exit(1);
    }
  }
}
