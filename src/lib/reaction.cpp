
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  *_gamma += *_beam - *_elec;
  // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)

  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);

  _emu_prime_mag2 = _elec->M2();
  _emu_mag2 = _beam->M2();
  _elec_energy = _elec->E();

  _elec_mom = _elec->P();
  _theta_e = _elec->Theta() * 180.0 / PI;

  if (_elec->Phi() > 0)
    _elec_phi = _elec->Phi() * 180 / PI;
  else if (_elec->Phi() < 0)
    _elec_phi = (_elec->Phi() + 2 * PI) * 180 / PI;
}
double Reaction::elec_prime_mass2() {
  if (_emu_prime_mag2 != _emu_prime_mag2) SetElec();

  return _emu_prime_mag2;
}
double Reaction::elec_mass2() {
  if (_emu_mag2 != _emu_mag2) SetElec();

  return _emu_mag2;
}

double Reaction::elec_E() {
  if (_elec_energy != _elec_energy) SetElec();
  // std::cout << " emec mom " << _elec_mom << std::endl;

  return _elec_energy;
}
double Reaction::elec_mom() {
  if (_elec_mom != _elec_mom) SetElec();
  // std::cout << " emec mom " << _elec_mom << std::endl;

  return _elec_mom;
}
double Reaction::elec_theta() {
  if (_theta_e != _theta_e) SetElec();

  return _theta_e;
}
double Reaction::elec_phi() {
  if (_elec_phi != _elec_phi) SetElec();

  return _elec_phi;
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}

void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->pid(i) == NEUTRON) {
    SetNeutron(i);
  } else {
    _numOther++;
    _hasOther = true;
    _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
}

float Reaction::w_hadron() {
  if (TwoPion_exclusive())
    return ((*_prot) + (*_pip) + (*_pim)).Mag();
  else
    return NAN;
}
float Reaction::w_difference() {
  if (TwoPion_exclusive())
    return (physics::W_calc(*_beam, *_elec) - ((*_prot) + (*_pip) + (*_pim)).Mag());
  else
    return NAN;
}

void Reaction::CalcMissMass() {
  auto mm_mpim = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();

  if (TwoPion_exclusive()) {
    *mm_excl += (*_gamma + *_target);
    *mm_excl -= *_prot;
    *mm_excl -= *_pip;
    *mm_excl -= *_pim;

    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();
    _excl_Mom = mm_excl->P();
  }
  if (TwoPion_missingPim()) {
    *mm_mpim += (*_gamma + *_target);

    *mm_mpim -= *_prot;
    *mm_mpim -= *_pip;
    _MM_mPim = mm_mpim->M();
    _MM2_mPim = mm_mpim->M2();
  }

  if (TwoPion_missingPip()) {
    *mm_mpip += (*_gamma + *_target);
    *mm_mpip -= *_prot;
    *mm_mpip -= *_pim;
    _MM2_mPip = mm_mpip->M2();
  }
  if (TwoPion_missingProt()) {
    *mm_mprot += (*_gamma + *_target);
    *mm_mprot -= *_pip;
    *mm_mprot -= *_pim;
    _MM2_mProt = mm_mprot->M2();
  }
}

float Reaction::MM2_exclusive() {
  if (_MM2_exclusive != _MM2_exclusive) CalcMissMass();
  return _MM2_exclusive;
}
float Reaction::MM_mPim() {
  if (_MM_mPim != _MM_mPim) CalcMissMass();
  return _MM_mPim;
}
float Reaction::MM2_mPim() {
  if (_MM2_mPim != _MM2_mPim) CalcMissMass();
  return _MM2_mPim;
}
float Reaction::MM2_mPip() {
  if (_MM2_mPip != _MM2_mPip) CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  if (_MM2_mProt != _MM2_mProt) CalcMissMass();
  return _MM2_mProt;
}
float Reaction::Energy_excl() {
  if (_excl_Energy != _excl_Energy) CalcMissMass();
  //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
  //  if (_x_mu_E > 0)
  return _excl_Energy;
  // else
  // return NAN;
}
float Reaction::Mom_excl() {
  if (_excl_Mom != _excl_Mom) CalcMissMass();
  return _excl_Mom;
  // else
  // return NAN;
}

MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_enrgy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_enrgy;
  _weight_mc = _data->mc_weight();
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  //_gamma = std::make_unique<TLorentzVector>();  // do i need this?
  _gamma_mc = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  //_elec = std::make_unique<TLorentzVector>();  // do i need this?
  _elec_mc = std::make_unique<TLorentzVector>();
  // this->SetElec();  // do i need this?
  this->SetMCElec();
  _prot_mc = std::make_unique<TLorentzVector>();
  _pip_mc = std::make_unique<TLorentzVector>();
  _pim_mc = std::make_unique<TLorentzVector>();
  //_other = std::make_unique<TLorentzVector>();  // do i need this?
  _other_mc = std::make_unique<TLorentzVector>();
  //_neutron = std::make_unique<TLorentzVector>();
}
// Reaction::~Reaction() {} // why this is not here
void MCReaction::SetMCElec() {
  //  _hasE = true;  //??
  _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);

  *_gamma_mc += *_beam - *_elec_mc;

  // Can calculate W and Q2 here
  _W_mc = physics::W_calc(*_beam, *_elec_mc);
  _Q2_mc = physics::Q2_calc(*_beam, *_elec_mc);
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
// }
float MCReaction::elec_mom_mc_gen() { return _elec_mc->P(); }
float MCReaction::elec_E_mc_gen() { return _elec_mc->E(); }

float MCReaction::elec_theta_mc_gen() { return _elec_mc->Theta() * 180 / PI; }
float MCReaction::elec_phi_mc_gen() {
  if (_elec_mc->Phi() >= 0)
    return (_elec_mc->Phi() * 180 / PI);
  else if (_elec_mc->Phi() < 0)
    return ((_elec_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}
