#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"

#include "podio/ROOTReader.h"
#include "podio/Frame.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"

#include "detectorSegmentations/GridDRcalo_k4geo.h"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  std::vector<std::string> filenames;

  for (int idx = 1; idx < argc; idx++)
    filenames.push_back( std::string(argv[idx]) );

  gStyle->SetOptFit(1);

  auto* pSeg = new dd4hep::DDSegmentation::GridDRcalo_k4geo("system:5,assembly:1,eta:-8,phi:9,x:32:-12,y:-12,c:1,module:2");

  TH1F* tN_C = new TH1F("N_C","Total Ceren p.e.;p.e.;Evt",400,0.,4000.);
  tN_C->Sumw2(); tN_C->SetLineColor(kBlue); tN_C->SetLineWidth(2);
  TH1F* tE_S = new TH1F("E_S","Energy of Scintillation ch.;GeV;Evt",1000,0.,5);
  tE_S->Sumw2(); tE_S->SetLineColor(kRed); tE_S->SetLineWidth(2);

  TH1F* tWav_C = new TH1F("wav_C","wavelength;nm;p.e.",120,300.,900.);
  tWav_C->Sumw2(); tWav_C->SetLineColor(kBlue); tWav_C->SetLineWidth(2);

  TH1F* tTime_C = new TH1F("time_C","time;ns;p.e.",150,10.,70.);
  tTime_C->Sumw2(); tTime_C->SetLineColor(kBlue); tTime_C->SetLineWidth(2);

  TH2F* tN2_C = new TH2F("N2_C",";x;y",150,-50.,100.,150,-50.,100.);
  TH2F* tN2_S = new TH2F("N2_S",";x;y",150,-50.,100.,150,-50.,100.);

  TH2F* tT2_C = new TH2F("T2_C",";nEta;nPhi",200,-100.,100.,200,0.,200.);
  TH2F* tT2_S = new TH2F("T2_S",";nEta;nPhi",200,-100.,100.,200,0.,200.);

  TH1F* tEnergy_C = new TH1F("energy_C","energy;GeV;Evt",1000,0.,100.);
  tEnergy_C->Sumw2(); tEnergy_C->SetLineColor(kBlue); tEnergy_C->SetLineWidth(2);   
  TH1F* tEnergy_S = new TH1F("energy_S","energy;GeV;Evt",1000,0.,100.);
  tEnergy_S->Sumw2(); tEnergy_S->SetLineColor(kRed); tEnergy_S->SetLineWidth(2);
  TH1F* tEnergy_Comb = new TH1F("energy_Comb","energy;GeV;Evt",2000,0.,200.);
  tEnergy_Comb->Sumw2(); tEnergy_Comb->SetLineColor(kGreen); tEnergy_Comb->SetLineWidth(2);

  double dCeren_Equalization_const = 0.016891891891892; // For new geometry
  double dScint_Equalization_const = 32.451728054518903; // For new geometry

//   double dCeren_SF = 1.0; // For new geometry
//   double dScint_SF = 1.0; // For new geometry

  auto reader = podio::ROOTReader();
  reader.openFiles(filenames);

  unsigned int entries = reader.getEntries("events");

  for (unsigned int iEvt = 0; iEvt < entries; iEvt++) {
    if (iEvt % 100 == 0) printf("Analyzing %dth event ...\n", iEvt);

    auto aframe = podio::Frame(reader.readEntry("events",iEvt));

    const auto* scintHits = static_cast<const edm4hep::SimCalorimeterHitCollection*>(aframe.get("DRcaloSiPMreadout_scint"));
    const auto* cerenHits = static_cast<const edm4hep::SimCalorimeterHitCollection*>(aframe.get("DRcaloSiPMreadoutSimHit"));

    const auto* rawTimeStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("DRcaloSiPMreadoutTimeStruct"));
    const auto* rawWavlenStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("DRcaloSiPMreadoutWaveLen"));

    //const auto* rawTimeStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("RawTimeStructs"));
    //const auto* rawWavlenStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("RawWavlenStructs"));
    //const auto* digiWaveforms = static_cast<const edm4hep::TimeSeriesCollection*>(aframe.get("DigiWaveforms"));
    //const auto* digiHits = static_cast<const edm4hep::CalorimeterHitCollection*>(aframe.get("DigiCalorimeterHits"));
    //const auto* procTimes = static_cast<const edm4hep::TimeSeriesCollection*>(aframe.get("DRpostprocTime"));
    //const auto* leakages = static_cast<const edm4hep::MCParticleCollection*>(aframe.get("Leakages"));

    double scintEn = 0.;
    int nCeren = 0;

    double dCeren_Calibrated_energy = 0.;
    double dScint_Calibrated_energy = 0.;

    /*for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
      const auto& scintHit = scintHits->at(idx);
      const auto cID = scintHit->getCellID();
      float en = scintHit.getEnergy();

      tT2_S->Fill(pSeg->numEta(cID),pSeg->numPhi(cID),en);
    }

    auto getMaxXY = [] (TH2* ahist) -> std::pair<int,int> {
      int x = 0;
      int y = 0;
      double val = 0.;

      for (unsigned ix=0; ahist->GetNbinsX()+2; ix++) {
        for (unsigned iy=0; ahist->GetNbinsY()+2; iy++) {
          if (ahist->GetBinContent(x,y) > val) {
            val = ahist->GetBinContent(x,y);
            x = ix;
            y = iy;
          }
        }
      }

      return std::make_pair(x,y);
    };

    auto*/

    for (unsigned int idx = 0; idx < scintHits->size(); idx++) {
      const auto& scintHit = scintHits->at(idx);
      const auto cID = scintHit.getCellID();
      float en = scintHit.getEnergy();

      tT2_S->Fill(pSeg->numEta(cID),pSeg->numPhi(cID),en);
      tN2_S->Fill(pSeg->x(cID),pSeg->y(cID),en);

      scintEn += en;
      dScint_Calibrated_energy += en * dScint_Equalization_const;
    }

    tE_S->Fill(scintEn);
    tEnergy_S->Fill(dScint_Calibrated_energy);
    for (unsigned int idx = 0; idx < cerenHits->size(); idx++) {
      const auto& cerenHit = cerenHits->at(idx);
      const auto cID = cerenHit.getCellID();
      int amp = cerenHit.getEnergy();

      tT2_C->Fill(pSeg->numEta(cID),pSeg->numPhi(cID),amp);
      tN2_C->Fill(pSeg->x(cID),pSeg->y(cID),amp);

      const auto& timeStruct = rawTimeStructs->at(idx);
      const auto& wavlenStruct = rawWavlenStructs->at(idx);

      for (unsigned int bin = 0; bin < timeStruct.adcCounts_size(); bin++) {
        int con = timeStruct.getAdcCounts(bin);
        tTime_C->Fill(timeStruct.getTime()+(static_cast<float>(bin)+0.5)*timeStruct.getInterval(),static_cast<float>(con));
      }

      for (unsigned int bin = 0; bin < wavlenStruct.adcCounts_size(); bin++) {
        int con = wavlenStruct.getAdcCounts(bin);
        tWav_C->Fill(wavlenStruct.getTime()+(static_cast<float>(bin)+0.5)*wavlenStruct.getInterval(),static_cast<float>(con));
      }

      nCeren += amp;
      dCeren_Calibrated_energy += amp * dCeren_Equalization_const;
    }

    tN_C->Fill(static_cast<float>(nCeren));
    tEnergy_C->Fill(dCeren_Calibrated_energy);

    tEnergy_Comb->Fill(dCeren_Calibrated_energy+dScint_Calibrated_energy);

    std::cout << "Evt: " << iEvt << " scintEn: " << scintEn << " nCeren: " << nCeren << " S ch E: " << dScint_Calibrated_energy << " C ch E: " << dCeren_Calibrated_energy << " Comb ch E: " << dScint_Calibrated_energy+dCeren_Calibrated_energy << std::endl;
  } // event loop

  std::string filenameStd = static_cast<std::string>(filenames.front());
  std::string extension = ".root";

  auto where = filenameStd.find(extension);
  if (where != std::string::npos) {
    filenameStd.erase(where, extension.length());
  }

  TString filename = static_cast<TString>(filenameStd);

  TFile* validFile = new TFile(filename+"_validation.root","RECREATE");

  validFile->WriteTObject(tN_C);
  validFile->WriteTObject(tE_S);
  validFile->WriteTObject(tEnergy_S);
  validFile->WriteTObject(tEnergy_C);
  validFile->WriteTObject(tEnergy_Comb);
  validFile->WriteTObject(tTime_C);
  validFile->WriteTObject(tWav_C);
  validFile->WriteTObject(tT2_S);
  validFile->WriteTObject(tN2_S);
  validFile->WriteTObject(tT2_C);
  validFile->WriteTObject(tN2_C);

  validFile->Close();
}
