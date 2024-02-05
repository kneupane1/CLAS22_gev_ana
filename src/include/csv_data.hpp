#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data
{

  float w;
  float w_mc;

  float q2_mc;
  float w_had;
  float w_diff;
  float q2;
  float gen_elec_E;
  float gen_elec_mom;
  float gen_elec_theta;
  float gen_elec_phi;

  float dc_R1_edge;
  float dc_R2_edge;
  float dc_R3_edge;

  float pip_theta_rec;
  float pip_phi_rec;
  int sec_pip;

  float weight_gen;

  // Static functions can be called without making a new struct
  static std::string header()
  {
    return "w_rec,q2_rec,dc_r1_edge,dc_r2_edge,dc_r3_edge,pip_theta_rec,pip_phi_rec,sec_pip,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data)
  {
    os << std::setprecision(7);
    // os << data.w_mc << ",";
    // os << data.q2_mc << ",";

    // // os << data.gen_elec_E << ",";
    // os << data.gen_elec_mom << ",";
    // os << data.gen_elec_theta << ",";
    // os << data.gen_elec_phi << ",";

    os << data.w << ",";
    os << data.q2 << ",";
    // os << data.w_had << ",";
    os << data.dc_R1_edge << ",";
    os << data.dc_R2_edge << ",";
    os << data.dc_R3_edge << ",";

    os << data.pip_theta_rec << ",";
    os << data.pip_phi_rec << ",";

    os << std::setprecision(1);
    os << data.sec_pip << ",";

    os << std::setprecision(7);
    os << data.weight_gen << ",";

    return os;
  };
};

#endif
