
#include "cuts.hpp"
#include <iostream>
#include "TFile.h"
#include "reaction.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12> &data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

bool Cuts::ElectronCuts()
{
  bool _elec = true;
  // Number of good particles is greater than 0
  // So that we can check at(0) without errors
  _elec &= (_data->gpart() > 0);
  if (!_elec)
    return false;

  _elec &= (_data->gpart() < 20);

  _elec &= (_data->charge(0) == NEGATIVE);
  _elec &= (_data->pid(0) == ELECTRON);
  _elec &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  return _elec;
}

bool Cuts::IsPip(int i)
{
  if (_data->gpart() <= i)
    return false;
  bool _pip = true;
  _pip &= (_data->charge(i) == POSITIVE);
  _pip &= (_data->pid(i) == PIP);
  _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
  // _pip &= ( abs(_dt->dt_ctof_Pi(i)) < 0.4);

  // _pip &= (4000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

  return _pip;
}
bool Cuts::IsProton(int i)
{
  if (_data->gpart() <= i)
    return false;
  bool _proton = true;
  _proton &= (_data->charge(i) == POSITIVE);
  _proton &= (_data->pid(i) == PROTON);
  _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
  // if (!std::isnan(_dt->dt_P(i))) _proton &= (abs(_dt->dt_ctof_P(i)) < 0.4);
  // _proton &= abs(_dt->dt_P(i)) < 0.5;
  // _proton &= abs(_dt->dt_ctof_P(i)) < 0.4;

  // _proton &= (4000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

  return _proton;
}
bool Cuts::IsPim(int i)
{
  if (_data->gpart() <= i)
    return false;
  bool _pim = true;
  _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (_data->pid(i) == PIM);
  _pim &= (abs(_dt->dt_Pi(i)) < 0.5) || (abs(_dt->dt_ctof_Pi(i)) < 0.5);
  // _pim &= (abs(_dt->dt_ctof_Pi(i)) < 0.5);
  // _pim &= (4000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

  return _pim;
}

bool rga_Cuts::ElectronCuts()
{
  bool cut = true;
  cut &= (_data->gpart() > 0);
  if (!cut)
    return false;

  cut &= (_data->gpart() < 20);
  //
  cut &= (_data->charge(0) == NEGATIVE);
  cut &= (_data->pid(0) == ELECTRON);
  cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);

  return cut;
}