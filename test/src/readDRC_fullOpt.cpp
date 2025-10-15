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

  TH1F* tN_C = new TH1F("N_C","Total Ceren p.e.;p.e.;Evt",1000,0.,10000.);
  tN_C->Sumw2(); tN_C->SetLineColor(kBlue); tN_C->SetLineWidth(2);
  TH1F* tN_S = new TH1F("N_S","Total Scintillation p.e.;p.e.;Evt",3000,0.,30000.);
  tN_S->Sumw2(); tN_S->SetLineColor(kRed); tN_S->SetLineWidth(2);

  TH1F* tWav_C = new TH1F("wav_C","wavelength;nm;p.e.",120,300.,900.);
  tWav_C->Sumw2(); tWav_C->SetLineColor(kBlue); tWav_C->SetLineWidth(2);
  TH1F* tTime_C = new TH1F("time_C","time;ns;p.e.",180,-10.,80.);
  tTime_C->Sumw2(); tTime_C->SetLineColor(kBlue); tTime_C->SetLineWidth(2);

  TH1F* tWav_S = new TH1F("wav_S","wavelength;nm;p.e.",120,300.,900.);
  tWav_S->Sumw2(); tWav_S->SetLineColor(kRed); tWav_S->SetLineWidth(2);
  TH1F* tTime_S = new TH1F("time_S","time;ns;p.e.",180,-10.,80.);
  tTime_S->Sumw2(); tTime_S->SetLineColor(kRed); tTime_S->SetLineWidth(2);

  TH2F* tN2_C = new TH2F("N2_C",";x;y",150,-50.,100.,150,-50.,100.);
  TH2F* tN2_S = new TH2F("N2_S",";x;y",150,-50.,100.,150,-50.,100.);

  TH2F* tT2_C = new TH2F("T2_C",";nEta;nPhi",80,-40.,40.,40,0.,40.);
  TH2F* tT2_S = new TH2F("T2_S",";nEta;nPhi",80,-40.,40.,40,0.,40.);

  TH1F* tEnergy_C = new TH1F("energy_C","energy;GeV;Evt",2000,0.,200.);
  tEnergy_C->Sumw2(); tEnergy_C->SetLineColor(kBlue); tEnergy_C->SetLineWidth(2);   
  TH1F* tEnergy_S = new TH1F("energy_S","energy;GeV;Evt",2000,0.,200.);
  tEnergy_S->Sumw2(); tEnergy_S->SetLineColor(kRed); tEnergy_S->SetLineWidth(2);
  TH1F* tEnergy_Comb = new TH1F("energy_Comb","energy;GeV;Evt",4000,0.,400.);
  tEnergy_Comb->Sumw2(); tEnergy_Comb->SetLineColor(kGreen); tEnergy_Comb->SetLineWidth(2);

  double dCeren_Equalization_const = 0.021376656690894; // For INCREADIBLE (barrel only) geometry
  double dScint_Equalization_const = 0.006568144499179; // For INCREADIBLE (barrel only) geometry

  auto reader = podio::ROOTReader();
  reader.openFiles(filenames);

  unsigned int entries = reader.getEntries("events");

  for (unsigned int iEvt = 0; iEvt < entries; iEvt++) {
    if (iEvt % 100 == 0) printf("Analyzing %dth event ...\n", iEvt);

    auto aframe = podio::Frame(reader.readEntry("events",iEvt));

    const auto* simHits = static_cast<const edm4hep::SimCalorimeterHitCollection*>(aframe.get("DRcaloSiPMreadoutSimHit"));
    const auto* rawTimeStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("DRcaloSiPMreadoutTimeStruct"));
    const auto* rawWavlenStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("DRcaloSiPMreadoutWaveLen"));

    //const auto* rawTimeStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("RawTimeStructs"));
    //const auto* rawWavlenStructs = static_cast<const edm4hep::RawTimeSeriesCollection*>(aframe.get("RawWavlenStructs"));
    //const auto* digiWaveforms = static_cast<const edm4hep::TimeSeriesCollection*>(aframe.get("DigiWaveforms"));
    //const auto* digiHits = static_cast<const edm4hep::CalorimeterHitCollection*>(aframe.get("DigiCalorimeterHits"));
    //const auto* procTimes = static_cast<const edm4hep::TimeSeriesCollection*>(aframe.get("DRpostprocTime"));
    //const auto* leakages = static_cast<const edm4hep::MCParticleCollection*>(aframe.get("Leakages"));

    int nScint = 0.;
    int nCeren = 0.;

    double dCeren_Calibrated_energy = 0.;
    double dScint_Calibrated_energy = 0.;
    double dComb_Calibrated_energy  = 0.;

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
    for (unsigned int idx = 0; idx < simHits->size(); idx++) {
      const auto& simHit = simHits->at(idx);
      const auto cID = simHit.getCellID();
      int amp = simHit.getEnergy();
      bool isCeren = pSeg->IsCerenkov(cID);

      if (isCeren) {
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
      else {
        tT2_S->Fill(pSeg->numEta(cID),pSeg->numPhi(cID),amp);
        tN2_S->Fill(pSeg->x(cID),pSeg->y(cID),amp);

        const auto& timeStruct = rawTimeStructs->at(idx);
        const auto& wavlenStruct = rawWavlenStructs->at(idx);

        for (unsigned int bin = 0; bin < timeStruct.adcCounts_size(); bin++) {
          int con = timeStruct.getAdcCounts(bin);
          tTime_S->Fill(timeStruct.getTime()+(static_cast<float>(bin)+0.5)*timeStruct.getInterval(),static_cast<float>(con));
        }

        for (unsigned int bin = 0; bin < wavlenStruct.adcCounts_size(); bin++) {
          int con = wavlenStruct.getAdcCounts(bin);
          tWav_S->Fill(wavlenStruct.getTime()+(static_cast<float>(bin)+0.5)*wavlenStruct.getInterval(),static_cast<float>(con));
        }

        nScint += amp;
        dScint_Calibrated_energy += amp * dScint_Equalization_const;
      }
    }

    dComb_Calibrated_energy = dScint_Calibrated_energy + dCeren_Calibrated_energy;

    std::cout << "Evt: " << iEvt << " nScint: " << nScint << " nCeren: " << nCeren << " S ch E: " << dScint_Calibrated_energy << " C ch E: " << dCeren_Calibrated_energy << std::endl;
    tN_S->Fill(static_cast<float>(nScint));
    tEnergy_S->Fill(dScint_Calibrated_energy);

    tN_C->Fill(static_cast<float>(nCeren));
    tEnergy_C->Fill(dCeren_Calibrated_energy);

    tEnergy_Comb->Fill(dComb_Calibrated_energy);
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
  validFile->WriteTObject(tN_S);
  validFile->WriteTObject(tEnergy_S);
  validFile->WriteTObject(tEnergy_C);
  validFile->WriteTObject(tEnergy_Comb);
  validFile->WriteTObject(tTime_C);
  validFile->WriteTObject(tWav_C);
  validFile->WriteTObject(tTime_S);
  validFile->WriteTObject(tWav_S);
  validFile->WriteTObject(tT2_S);
  validFile->WriteTObject(tN2_S);
  validFile->WriteTObject(tT2_C);
  validFile->WriteTObject(tN2_C);

  validFile->Close();
}
