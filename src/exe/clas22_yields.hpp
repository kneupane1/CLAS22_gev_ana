
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  float beam_energy = 10.6;
  if (std::is_same<CutType, rga_Cuts>::value) {
    beam_energy = 10.6;
  } else if (std::is_same<CutType, uconn_Cuts>::value) {
    beam_energy = 10.6;
  }

  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  // for sim data use it
  auto data = std::make_shared<Branches12>(_chain, true);
  // for exp data use it
  // auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  float vertex_hadron[3][3];

  size_t total_twopion_events = 0;

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    int statusPim = -9999;
    int statusPip = -9999;
    int statusProt = -9999;

    if (data->mc_npart() < 1) continue;

    // // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

    for (int part = 1; part < data->mc_npart(); part++) {
      // Check particle ID's and fill the reaction class

      if (data->mc_pid(part) == PIP) {
        mc_event->SetMCPip(part);

      } else if (data->mc_pid(part) == PROTON) {
        mc_event->SetMCProton(part);

      } else if (data->mc_pid(part) == PIM) {
        mc_event->SetMCPim(part);

        // } else {
        //   mc_event->SetMCOther(part);
      }
    }

    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<uconn_Cuts>(data);
    // auto cuts = std::make_shared<rga_Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);

    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);

      // Check particle ID's and fill the reaction class
      if (cuts->IsProton(part)) {
        event->SetProton(part);
        statusProt = abs(data->status(part));

      } else if (cuts->IsPip(part)) {
        event->SetPip(part);
        statusPip = abs(data->status(part));

      } else if (cuts->IsPim(part)) {
        event->SetPim(part);
        statusPim = abs(data->status(part));

      } else {
        event->SetOther(part);
      }
    }

    if (event->TwoPion_exclusive()) {
      if (event->W() > 1.25 && event->W() < 2.55 && event->Q2() > 1.5 && event->Q2() < 30.5) {

        total_twopion_events++;
        csv_data output;

        // // // // // //  1) for generated
        output.w_mc = mc_event->W_mc();
        output.q2_mc = mc_event->Q2_mc();

        // // // /// 2) reconstructed  and rec exclusive
        output.w = event->W();
        output.q2 = event->Q2();
        output.w_had = event->w_hadron();

        output.weight_gen = event->weight();

        _sync->write(output);
      }
    }
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  std::cout << " total no of events = " << total << std::endl;

  return num_of_events;
}
#endif
