// random.h

#ifndef RANDOM_H
#define RANDOM_H

typedef struct Random {
  int seed;
  void* randomObject;
} Random;

Random* createRandom(unsigned long seed); // Set up a random generator
void deleteRandom(Random* rand); // Delete the random generator
double randDouble(Random* rand); // Return a random float in [0,1)
int randInt(Random*, int min, int max);// Return a random int in [min, max)

#endif
