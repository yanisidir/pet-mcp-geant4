#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

void analyze_coincidence_events(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* eventTree = nullptr;
  file->GetObject("EventSummaryTree", eventTree);
  if (!eventTree) {
    std::cerr << "Missing EventSummaryTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t plusCount = 0;
  Int_t minusCount = 0;
  Bool_t isCoincidence = false;

  eventTree->SetBranchAddress("electronChannelPlusCount",
                              &plusCount);
  eventTree->SetBranchAddress("electronChannelMinusCount",
                              &minusCount);
  eventTree->SetBranchAddress("isCoincidence",
                              &isCoincidence);

  const Long64_t totalEventCount = eventTree->GetEntries();
  std::vector<Int_t> plusMultiplicities;
  std::vector<Int_t> minusMultiplicities;

  Int_t maxPlusCount = 0;
  Int_t maxMinusCount = 0;

  for (Long64_t entry = 0; entry < totalEventCount; ++entry) {
    eventTree->GetEntry(entry);
    if (!isCoincidence) {
      continue;
    }

    plusMultiplicities.push_back(plusCount);
    minusMultiplicities.push_back(minusCount);
    maxPlusCount = std::max(maxPlusCount, plusCount);
    maxMinusCount = std::max(maxMinusCount, minusCount);
  }

  const Long64_t coincidenceEventCount =
    static_cast<Long64_t>(plusMultiplicities.size());

  std::cout << "File: " << fileName << std::endl;
  if (coincidenceEventCount == 0) {
    std::cout << "\nNombre total d'evenements : "
              << totalEventCount << std::endl;
    std::cout << "Nombre d'evenements en coincidence : 0"
              << std::endl;
    std::cout << "No coincidence event available for analysis."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  const Int_t histogramMaximum =
    std::max(maxPlusCount, maxMinusCount);
  const Int_t histogramBinCount = histogramMaximum + 1;

  TH1D* plusHistogram =
    new TH1D("hCoincidencePlusMultiplicity",
             "Electron multiplicity in coincidence events;"
             "Electrons reaching MCP channels;Events",
             histogramBinCount,
             -0.5,
             histogramMaximum + 0.5);

  TH1D* minusHistogram =
    new TH1D("hCoincidenceMinusMultiplicity",
             "Electron multiplicity in coincidence events;"
             "Electrons reaching MCP channels;Events",
             histogramBinCount,
             -0.5,
             histogramMaximum + 0.5);

  TH2D* correlationHistogram =
    new TH2D("hCoincidencePlusVsMinus",
             "Electron multiplicity correlation in coincidence events;"
             "N(+z);N(-z)",
             maxPlusCount + 1,
             -0.5,
             maxPlusCount + 0.5,
             maxMinusCount + 1,
             -0.5,
             maxMinusCount + 0.5);

  Double_t sumPlus = 0.0;
  Double_t sumMinus = 0.0;
  Double_t sumPlusSquared = 0.0;
  Double_t sumMinusSquared = 0.0;
  Double_t sumPlusMinus = 0.0;
  Double_t sumAbsoluteDifference = 0.0;

  for (std::size_t eventIndex = 0;
       eventIndex < plusMultiplicities.size();
       ++eventIndex) {
    const Double_t plus =
      static_cast<Double_t>(plusMultiplicities[eventIndex]);
    const Double_t minus =
      static_cast<Double_t>(minusMultiplicities[eventIndex]);

    plusHistogram->Fill(plus);
    minusHistogram->Fill(minus);
    correlationHistogram->Fill(plus, minus);

    sumPlus += plus;
    sumMinus += minus;
    sumPlusSquared += plus*plus;
    sumMinusSquared += minus*minus;
    sumPlusMinus += plus*minus;
    sumAbsoluteDifference += std::fabs(plus - minus);
  }

  const Double_t eventCount =
    static_cast<Double_t>(coincidenceEventCount);
  const Double_t meanPlus = sumPlus/eventCount;
  const Double_t meanMinus = sumMinus/eventCount;
  const Double_t meanAbsoluteDifference =
    sumAbsoluteDifference/eventCount;

  const Double_t variancePlus =
    sumPlusSquared/eventCount - meanPlus*meanPlus;
  const Double_t varianceMinus =
    sumMinusSquared/eventCount - meanMinus*meanMinus;
  const Double_t covariance =
    sumPlusMinus/eventCount - meanPlus*meanMinus;

  Double_t correlationCoefficient =
    std::numeric_limits<Double_t>::quiet_NaN();
  if (variancePlus > 0.0 && varianceMinus > 0.0) {
    correlationCoefficient =
      covariance/std::sqrt(variancePlus*varianceMinus);
  }

  gStyle->SetOptStat(0);

  plusHistogram->SetLineColor(kBlue + 1);
  plusHistogram->SetLineWidth(3);
  plusHistogram->SetMarkerColor(kBlue + 1);
  plusHistogram->SetMarkerStyle(20);
  plusHistogram->SetMarkerSize(0.8);

  minusHistogram->SetLineColor(kGreen + 2);
  minusHistogram->SetLineWidth(3);
  minusHistogram->SetMarkerColor(kGreen + 2);
  minusHistogram->SetMarkerStyle(21);
  minusHistogram->SetMarkerSize(0.8);

  const Double_t displayMaximum =
    1.15*std::max(plusHistogram->GetMaximum(),
                  minusHistogram->GetMaximum());
  plusHistogram->SetMaximum(displayMaximum);

  TCanvas* multiplicityCanvas =
    new TCanvas("cCoincidenceMultiplicity",
                "Coincidence electron multiplicity",
                950,
                720);
  multiplicityCanvas->SetGridy();
  plusHistogram->Draw("HIST E");
  minusHistogram->Draw("HIST E SAME");

  TLegend* legend = new TLegend(0.66, 0.73, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(plusHistogram, "+z MCP stack", "l");
  legend->AddEntry(minusHistogram, "-z MCP stack", "l");
  legend->Draw();

  multiplicityCanvas->SaveAs(
    "coincidence_electron_multiplicity.png");

  gStyle->SetPaintTextFormat("g");
  TCanvas* correlationCanvas =
    new TCanvas("cCoincidencePlusVsMinus",
                "Coincidence plus versus minus",
                900,
                760);
  correlationCanvas->SetRightMargin(0.15);
  correlationCanvas->SetGrid();
  correlationHistogram->Draw("COLZ");
  correlationCanvas->SaveAs("coincidence_plus_vs_minus.png");

  TH2D* percentHistogram =
    static_cast<TH2D*>(
      correlationHistogram->Clone("hCoincidencePlusVsMinusPercent"));
  percentHistogram->SetTitle(
    "Electron multiplicity correlation in coincidence events (%);"
    "N(+z);N(-z)");
  percentHistogram->Scale(100.0/eventCount);

  gStyle->SetPaintTextFormat("4.1f");
  TCanvas* percentCanvas =
    new TCanvas("cCoincidencePlusVsMinusPercent",
                "Coincidence plus versus minus percent",
                900,
                760);
  percentCanvas->SetRightMargin(0.15);
  percentCanvas->SetGrid();
  percentHistogram->SetMarkerSize(1.3);
  percentHistogram->Draw("COLZ TEXT");
  percentCanvas->SaveAs(
    "coincidence_plus_vs_minus_percent.png");

  const Int_t plusOneBin =
    correlationHistogram->GetXaxis()->FindBin(1.0);
  const Int_t minusOneBin =
    correlationHistogram->GetYaxis()->FindBin(1.0);
  const Double_t oneOneCount =
    correlationHistogram->GetBinContent(plusOneBin, minusOneBin);
  const Double_t oneOneFraction = oneOneCount/eventCount;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nStatistiques des evenements en coincidence"
            << std::endl;
  std::cout << "Nombre total d'evenements : "
            << totalEventCount << std::endl;
  std::cout << "Nombre d'evenements en coincidence : "
            << coincidenceEventCount << std::endl;

  std::cout << "\nMean N(+z) : " << meanPlus << std::endl;
  std::cout << "RMS N(+z) : "
            << plusHistogram->GetRMS() << std::endl;
  std::cout << "Max N(+z) : " << maxPlusCount << std::endl;

  std::cout << "\nMean N(-z) : " << meanMinus << std::endl;
  std::cout << "RMS N(-z) : "
            << minusHistogram->GetRMS() << std::endl;
  std::cout << "Max N(-z) : " << maxMinusCount << std::endl;

  std::cout << "\nMean |N(+z)-N(-z)| : "
            << meanAbsoluteDifference << std::endl;
  std::cout << "Correlation coefficient : ";
  if (std::isnan(correlationCoefficient)) {
    std::cout << "undefined (zero variance)";
  } else {
    std::cout << correlationCoefficient;
  }
  std::cout << std::endl;

  std::cout << "TH2 correlation factor : "
            << correlationHistogram->GetCorrelationFactor()
            << std::endl;
  std::cout << "Fraction with N(+z)=1 and N(-z)=1 : "
            << oneOneFraction
            << " (" << 100.0*oneOneFraction << " %, "
            << static_cast<Long64_t>(oneOneCount)
            << " events)" << std::endl;

  std::cout << "\nSaved: coincidence_electron_multiplicity.png"
            << std::endl;
  std::cout << "Saved: coincidence_plus_vs_minus.png"
            << std::endl;
  std::cout << "Saved: coincidence_plus_vs_minus_percent.png"
            << std::endl;

  file->Close();
  delete file;
}
