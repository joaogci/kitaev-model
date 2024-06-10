#include "replica.h"

void init_replica(Flux* flux, double beta, Replica *replica)
{
  replica->flux = flux;
  replica->beta = beta;
}

void exchange(double r, Replica *replica_i, Replica *replica_j)
{
  double log_f;
  Flux *flux_tmp;

  log_f = log_weight(replica_i->beta, replica_j->flux) + log_weight(replica_j->beta, replica_i->flux);
  log_f += - log_weight(replica_i->beta, replica_i->flux) - log_weight(replica_j->beta, replica_j->flux);

  if (r <= exp(log_f)) {
    flux_tmp = replica_i->flux;
    replica_i->flux = replica_j->flux;
    replica_j->flux = flux_tmp;
  }
}
