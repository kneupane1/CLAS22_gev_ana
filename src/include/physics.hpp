
#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include <TLorentzVector.h>
#include "TROOT.h"
#include "constants.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime);
}  // namespace physics

#endif
