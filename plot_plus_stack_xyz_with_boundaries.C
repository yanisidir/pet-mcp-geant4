#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{
const Int_t kPlusSide = 1;
const Int_t kNumberOfMCPs = 3;

const Double_t kMcpLengthCm = 0.3;
const Double_t kMcpGapCm = 0.1;
const Double_t kFirstMcpCenterZCm = 20.0;

const Double_t kXMinimumCm = -0.3;
const Double_t kXMaximumCm = 0.3;
const Double_t kYMinimumCm = -0.3;
const Double_t kYMaximumCm = 0.3;
const Double_t kZMinimumCm = 19.8;
const Double_t kZMaximumCm = 21.0;

struct PlateGraph
{
  Int_t mcpIndex;
  Color_t color;
  TGraph2D* graph;
  Long64_t totalEntries;
};

std::string FindBranch(TTree* tree,
                       const std::vector<std::string>& candidates)
{
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    if (tree->GetBranch(candidates[i].c_str())) {
      return candidates[i];
    }
  }
  return "";
}

TLeaf* FindLeaf(TTree* tree,
                const std::vector<std::string>& candidates,
                std::string& selectedName)
{
  selectedName = FindBranch(tree, candidates);
  if (selectedName.empty()) {
    return nullptr;
  }
  return tree->GetLeaf(selectedName.c_str());
}

PlateGraph* FindPlate(std::vector<PlateGraph>& plates, Int_t mcpIndex)
{
  for (std::size_t i = 0; i < plates.size(); ++i) {
    if (plates[i].mcpIndex == mcpIndex) {
      return &plates[i];
    }
  }
  return nullptr;
}

void AddLine(std::vector<TPolyLine3D*>& lines,
             Double_t x1,
             Double_t y1,
             Double_t z1,
             Double_t x2,
             Double_t y2,
             Double_t z2,
             Color_t color,
             Style_t style = 2)
{
  TPolyLine3D* line = new TPolyLine3D(2);
  line->SetPoint(0, x1, y1, z1);
  line->SetPoint(1, x2, y2, z2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  line->Draw();
  lines.push_back(line);
}

void DrawMcpWireframe(std::vector<TPolyLine3D*>& lines,
                      Double_t zMinimum,
                      Double_t zMaximum,
                      Color_t color)
{
  const Double_t x[4] = {
    kXMinimumCm, kXMaximumCm, kXMaximumCm, kXMinimumCm
  };
  const Double_t y[4] = {
    kYMinimumCm, kYMinimumCm, kYMaximumCm, kYMaximumCm
  };

  for (Int_t corner = 0; corner < 4; ++corner) {
    const Int_t nextCorner = (corner + 1) % 4;

    AddLine(lines,
            x[corner], y[corner], zMinimum,
            x[nextCorner], y[nextCorner], zMinimum,
            color);
    AddLine(lines,
            x[corner], y[corner], zMaximum,
            x[nextCorner], y[nextCorner], zMaximum,
            color);
    AddLine(lines,
            x[corner], y[corner], zMinimum,
            x[corner], y[corner], zMaximum,
            color);
  }
}
}  // namespace

void plot_plus_stack_xyz_with_boundaries(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* tree = nullptr;
  file->GetObject("ElectronChannelHitTree", tree);
  if (!tree) {
    std::cerr << "ElectronChannelHitTree was not found in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::string xName;
  std::string yName;
  std::string zName;
  std::string sideName;
  std::string mcpIndexName;

  TLeaf* xLeaf = FindLeaf(
    tree, {"x_cm", "x", "posX", "positionX"}, xName);
  TLeaf* yLeaf = FindLeaf(
    tree, {"y_cm", "y", "posY", "positionY"}, yName);
  TLeaf* zLeaf = FindLeaf(
    tree, {"z_cm", "z", "posZ", "positionZ"}, zName);
  TLeaf* sideLeaf = FindLeaf(tree, {"side"}, sideName);
  TLeaf* mcpIndexLeaf = FindLeaf(
    tree, {"mcpIndex", "McpIndex", "plateIndex"}, mcpIndexName);

  if (!xLeaf || !yLeaf || !zLeaf || !sideLeaf || !mcpIndexLeaf) {
    std::cerr << "Required branches could not be detected." << std::endl;
    std::cerr << "  position: " << xName << ", " << yName << ", "
              << zName << std::endl;
    std::cerr << "  plate: " << sideName << ", " << mcpIndexName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "ElectronChannelHitTree branch detection:" << std::endl;
  std::cout << "  position: " << xName << ", " << yName << ", "
            << zName << std::endl;
  std::cout << "  filter: " << sideName << " == 1" << std::endl;
  std::cout << "  plate index: " << mcpIndexName << std::endl;

  const Color_t colors[kNumberOfMCPs] = {
    kRed + 1, kOrange + 7, kBlue + 1
  };

  std::vector<PlateGraph> plates;
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    PlateGraph plate;
    plate.mcpIndex = mcpIndex;
    plate.color = colors[mcpIndex];
    plate.graph = new TGraph2D();
    plate.graph->SetName(
      ("plusStackMcp" + std::to_string(mcpIndex)).c_str());
    plate.totalEntries = 0;
    plates.push_back(plate);
  }

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);

    if (static_cast<Int_t>(sideLeaf->GetValue()) != kPlusSide) {
      continue;
    }

    PlateGraph* plate = FindPlate(
      plates, static_cast<Int_t>(mcpIndexLeaf->GetValue()));
    if (!plate) {
      continue;
    }
    plate->totalEntries++;

    const Double_t x = xLeaf->GetValue();
    const Double_t y = yLeaf->GetValue();
    const Double_t z = zLeaf->GetValue();

    if (x < kXMinimumCm || x > kXMaximumCm ||
        y < kYMinimumCm || y > kYMaximumCm ||
        z < kZMinimumCm || z > kZMaximumCm) {
      continue;
    }

    plate->graph->SetPoint(plate->graph->GetN(), x, y, z);
  }

  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas(
    "plusStackXyzCanvas",
    "Electron positions in the +z MCP stack",
    1200,
    900);
  canvas->SetTheta(22.0);
  canvas->SetPhi(38.0);
  canvas->SetLeftMargin(0.10);
  canvas->SetRightMargin(0.20);

  // Invisible corner points impose the requested presentation ranges.
  TGraph2D* frame = new TGraph2D();
  frame->SetName("plusStackXyzFrame");
  for (Int_t xIndex = 0; xIndex < 2; ++xIndex) {
    for (Int_t yIndex = 0; yIndex < 2; ++yIndex) {
      for (Int_t zIndex = 0; zIndex < 2; ++zIndex) {
        frame->SetPoint(
          frame->GetN(),
          xIndex == 0 ? kXMinimumCm : kXMaximumCm,
          yIndex == 0 ? kYMinimumCm : kYMaximumCm,
          zIndex == 0 ? kZMinimumCm : kZMaximumCm);
      }
    }
  }

  frame->SetTitle(
    "Electron positions in the +z MCP stack;x [cm];y [cm];z [cm]");
  frame->SetMarkerStyle(1);
  frame->SetMarkerSize(0.0);
  frame->Draw("P");

  TLegend* legend = new TLegend(0.14, 0.70, 0.37, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    plates[plate].graph->SetMarkerStyle(20);
    plates[plate].graph->SetMarkerSize(0.35);
    plates[plate].graph->SetMarkerColor(plates[plate].color);
    if (plates[plate].graph->GetN() > 0) {
      plates[plate].graph->Draw("P SAME");
    }

    const std::string label =
      "+z MCP " + std::to_string(plates[plate].mcpIndex);
    legend->AddEntry(plates[plate].graph, label.c_str(), "p");
  }

  std::vector<TPolyLine3D*> boundaryLines;
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const Double_t centerZ =
      kFirstMcpCenterZCm +
      mcpIndex*(kMcpLengthCm + kMcpGapCm);
    const Double_t minimumZ = centerZ - 0.5*kMcpLengthCm;
    const Double_t maximumZ = centerZ + 0.5*kMcpLengthCm;

    DrawMcpWireframe(
      boundaryLines, minimumZ, maximumZ, colors[mcpIndex]);
  }

  legend->Draw();

  TLatex labels;
  labels.SetNDC();
  labels.SetTextSize(0.027);
  labels.SetTextAlign(12);
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const Double_t centerZ =
      kFirstMcpCenterZCm +
      mcpIndex*(kMcpLengthCm + kMcpGapCm);
    const Double_t minimumZ = centerZ - 0.5*kMcpLengthCm;
    const Double_t maximumZ = centerZ + 0.5*kMcpLengthCm;

    labels.SetTextColor(colors[mcpIndex]);
    labels.DrawLatex(
      0.71,
      0.86 - 0.045*mcpIndex,
      Form("MCP%d: %.2f < z < %.2f cm",
           mcpIndex, minimumZ, maximumZ));
  }

  canvas->Modified();
  canvas->Update();
  canvas->SaveAs("plus_stack_xyz_with_mcp_boundaries.png");

  std::cout << "\nElectrons in side +z:" << std::endl;
  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    std::cout << "  MCP " << plates[plate].mcpIndex
              << ": total=" << plates[plate].totalEntries
              << ", displayed=" << plates[plate].graph->GetN()
              << std::endl;
  }
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "\nMCP z boundaries [cm]:" << std::endl;
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const Double_t centerZ =
      kFirstMcpCenterZCm +
      mcpIndex*(kMcpLengthCm + kMcpGapCm);
    std::cout << "  MCP " << mcpIndex << ": "
              << centerZ - 0.5*kMcpLengthCm << " -> "
              << centerZ + 0.5*kMcpLengthCm << std::endl;
  }
  std::cout << "\nSaved: plus_stack_xyz_with_mcp_boundaries.png"
            << std::endl;

  file->Close();
  delete file;
}
