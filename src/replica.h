#ifndef REPLICA_H
#define REPLICA_H

#include <stdio.h>
#include <stdlib.h>

#include "monte_carlo.h"
#include "flux.h"

typedef struct Replica {
  Flux *flux;
  double beta;
} Replica;

/*
 * Initializes Replica struct
 */
void init_replica(Flux* flux, double beta, Replica *replica);

/*
 * Exchanges two replicas of the system
 */
void exchange(double r, Replica *replica_i, Replica *replica_j);

#endif // REPLICA_H
