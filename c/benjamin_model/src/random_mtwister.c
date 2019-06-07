// random_mtwister.c
// Wrapper file to a meserenne twister random generator

#include <stdlib.h>
#include <time.h>
#include "random.h"
#include "mtwister.h"

Random* createRandom(unsigned long seed) {
  Random* random = malloc(sizeof *random);
  random->seed = seed;
  srand(seed);
  random->randomObject = malloc(sizeof(MTRand));
  *((MTRand*) random->randomObject) = seedRand(rand());
  return random;
}

void deleteRandom(Random* random) {
  free((MTRand*)random->randomObject);
  free(random);
}

double randDouble(Random* random) {
  return genRand((MTRand*)random->randomObject);
}

int randInt(Random* random, int min, int max) {
  return min + (int) ((max-min) * genRand((MTRand*)random->randomObject));
}
