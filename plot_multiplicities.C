#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>

void plot_multiplicities(const char* fileName = "build/mcp_output_t0.root")
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
    std::cerr << "Missing EventSummaryTree in " << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t electronProducedCount = 0;
  Int_t electronChannelCount = 0;

  eventTree->SetBranchAddress("electronProducedCount",
                              &electronProducedCount);
  eventTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);

  const Long64_t nEntries = eventTree->GetEntries();

  Int_t maxProduced = 0;
  Int_t maxChannel = 0;

  for (Long64_t i = 0; i < nEntries; ++i) {
    eventTree->GetEntry(i);

    if (electronProducedCount > maxProduced) {
      maxProduced = electronProducedCount;
    }

    if (electronChannelCount > maxChannel) {
      maxChannel = electronChannelCount;
    }
  }

  const Int_t maxMultiplicity =
    std::max(maxProduced, maxChannel);

  TH1I* hProduced =
    new TH1I("hProduced",
             "Electron multiplicity per event;Number of electrons;Number of events",
             maxMultiplicity + 1,
             -0.5,
             maxMultiplicity + 0.5);

  TH1I* hChannel =
    new TH1I("hChannel",
             "Channel electron multiplicity per event;Number of electrons;Number of events",
             maxMultiplicity + 1,
             -0.5,
             maxMultiplicity + 0.5);

  for (Long64_t i = 0; i < nEntries; ++i) {
    eventTree->GetEntry(i);
    hProduced->Fill(electronProducedCount);
    hChannel->Fill(electronChannelCount);
  }

  gStyle->SetOptStat(1110);

  TCanvas* c1 = new TCanvas("c1",
                            "Produced electron multiplicity",
                            900,
                            700);
  hProduced->SetLineWidth(2);
  hProduced->Draw("HIST");
  c1->SaveAs("produced_electron_multiplicity.png");

  TCanvas* c2 = new TCanvas("c2",
                            "Channel electron multiplicity",
                            900,
                            700);
  hChannel->SetLineWidth(2);
  hChannel->Draw("HIST");
  c2->SaveAs("channel_electron_multiplicity.png");

  TCanvas* c3 = new TCanvas("c3",
                            "Comparison of electron multiplicities",
                            900,
                            700);
  hProduced->SetLineWidth(2);
  hChannel->SetLineWidth(2);
  hChannel->SetLineColor(kRed);
  hProduced->Draw("HIST");
  hChannel->Draw("HIST SAME");

  TLegend* legend = new TLegend(0.60, 0.75, 0.88, 0.88);
  legend->AddEntry(hProduced, "Produced electrons", "l");
  legend->AddEntry(hChannel, "Electrons reaching MCP channel", "l");
  legend->Draw();

  c3->SaveAs("electron_multiplicity_comparison.png");

  file->Close();
  delete file;
}