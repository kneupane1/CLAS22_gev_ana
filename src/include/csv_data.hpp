#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {

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


  float weight_gen;


  // Static functions can be called without making a new struct
  static std::string header() {
    return "w_gen,q2_gen,w_rec,q2_rec,w_had,weight";

  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    os << std::setprecision(7);
    os << data.w_mc << ",";
    os << data.q2_mc << ",";

    // // os << data.gen_elec_E << ",";
    // os << data.gen_elec_mom << ",";
    // os << data.gen_elec_theta << ",";
    // os << data.gen_elec_phi << ",";

    os << data.w << ",";
    os << data.q2 << ",";
    os << data.w_had << ",";

    os << data.weight_gen<< ",";

    return os;
  };
};

#endif
