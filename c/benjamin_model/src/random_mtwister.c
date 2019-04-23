// random_mtwister.c
// Wrapper file to a meserenne twister random generator

#include <stdlib.h>
#include <time.h>
#include "random.h"
#include "mtwister.h"

Random* createRandom(long _seed) {
  Random* random = malloc(sizeof *random);
  random->seed = _seed;
  srand(_seed);
  MTRand* mtrand = malloc(sizeof *mtrand);
  *mtrand = seedRand(rand());
  random->randomObject = mtrand;
  return random;
}

void deleteRandom(Random* random) {
  free(random->randomObject);
  free(random);
}

double randDouble(Random* random) {
  return genRand((MTRand*)random->randomObject);
}

int randInt(Random* random, int min, int max) {
  return min + (int) ((max-min) * genRand((MTRand*)random->randomObject));
}
