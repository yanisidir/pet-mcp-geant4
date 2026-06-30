#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{
typedef std::pair<int, int> GammaKey;

TLeaf* RequireLeaf(TTree* tree, const char* name)
{
  TLeaf* leaf = tree ? tree->GetLeaf(name) : nullptr;
  if (!leaf) {
    std::cerr << "Missing branch: " << name << std::endl;
  }
  return leaf;
}

void LoadGammaTrackInfo(TTree* tree,
                        std::set<GammaKey>& gammaTrackKeys,
                        std::map<GammaKey, std::set<int> >& parentIDsByGamma,
                        std::map<GammaKey, std::vector<double> >& energiesByGamma)
{
  if (!tree) {
    return;
  }

  TLeaf* eventLeaf = tree->GetLeaf("eventID");
  TLeaf* trackLeaf = tree->GetLeaf("trackID");
  TLeaf* parentLeaf = tree->GetLeaf("parentID");
  TLeaf* energyLeaf = tree->GetLeaf("kineticEnergy_keV");

  if (!eventLeaf || !trackLeaf || !parentLeaf) {
    return;
  }

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    const GammaKey key(
      static_cast<int>(eventLeaf->GetValue()),
      static_cast<int>(trackLeaf->GetValue()));
    const int parentID = static_cast<int>(parentLeaf->GetValue());

    gammaTrackKeys.insert(key);
    parentIDsByGamma[key].insert(parentID);

    if (energyLeaf) {
      energiesByGamma[key].push_back(energyLeaf->GetValue());
    }
  }
}
}

void DiagnoseElectronCreationInfo(
  const char* fileName = "build/mcp_output.root")
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
    std::cerr << "Missing ElectronChannelHitTree in " << fileName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  TLeaf* hasValidLeaf = RequireLeaf(tree, "hasValidCreationInfo");
  TLeaf* primaryGammaLeaf = RequireLeaf(tree, "primaryGammaTrackID");
  TLeaf* parentGammaLeaf = RequireLeaf(tree, "parentGammaTrackID");
  TLeaf* processLeaf = RequireLeaf(tree, "creatorProcessName");
  TLeaf* eventLeaf = RequireLeaf(tree, "eventID");

  if (!hasValidLeaf || !primaryGammaLeaf ||
      !parentGammaLeaf || !processLeaf || !eventLeaf) {
    std::cerr << "\nThis ROOT file was probably produced before the "
              << "latest creation-info patch." << std::endl;
    tree->Print();
    file->Close();
    delete file;
    return;
  }

  char creatorProcessName[64];
  creatorProcessName[0] = '\0';
  tree->SetBranchAddress("creatorProcessName", creatorProcessName);

  Long64_t total = tree->GetEntries();
  Long64_t valid = 0;
  Long64_t invalid = 0;
  Long64_t validPrimaryUnknown = 0;
  Long64_t validParentUnknown = 0;
  std::map<std::string, Long64_t> invalidByProcess;
  std::map<std::string, Long64_t> validPrimaryUnknownByProcess;
  std::map<std::string, Long64_t> validByProcess;
  std::map<int, Long64_t> validPrimaryUnknownByParentGamma;
  std::map<GammaKey, Long64_t> unknownPrimaryRowsByGamma;
  std::set<GammaKey> unknownPrimaryGammaKeys;

  for (Long64_t entry = 0; entry < total; ++entry) {
    tree->GetEntry(entry);

    const bool hasValid =
      static_cast<int>(hasValidLeaf->GetValue()) != 0;
    const int primaryGammaTrackID =
      static_cast<int>(primaryGammaLeaf->GetValue());
    const int parentGammaTrackID =
      static_cast<int>(parentGammaLeaf->GetValue());
    const int eventID = static_cast<int>(eventLeaf->GetValue());
    const std::string processName = creatorProcessName;

    if (hasValid) {
      ++valid;
      ++validByProcess[processName];

      if (primaryGammaTrackID < 0) {
        ++validPrimaryUnknown;
        ++validPrimaryUnknownByProcess[processName];
        ++validPrimaryUnknownByParentGamma[parentGammaTrackID];
        ++unknownPrimaryRowsByGamma[
          GammaKey(eventID, parentGammaTrackID)];
        unknownPrimaryGammaKeys.insert(
          GammaKey(eventID, parentGammaTrackID));
      }
      if (parentGammaTrackID < 0) {
        ++validParentUnknown;
      }
    } else {
      ++invalid;
      ++invalidByProcess[processName];
    }
  }

  std::cout << "\n=== Electron creation-info diagnostic ==="
            << std::endl;
  std::cout << "File: " << fileName << std::endl;
  std::cout << "Total ElectronChannelHitTree rows: " << total
            << std::endl;
  std::cout << "Rows with hasValidCreationInfo == true : " << valid
            << std::endl;
  std::cout << "Rows with hasValidCreationInfo == false: " << invalid
            << std::endl;
  std::cout << "Valid rows with primaryGammaTrackID == -1: "
            << validPrimaryUnknown << std::endl;
  std::cout << "Valid rows with parentGammaTrackID == -1 : "
            << validParentUnknown << std::endl;

  std::cout << "\ncreatorProcessName for invalid rows" << std::endl;
  if (invalidByProcess.empty()) {
    std::cout << "  none" << std::endl;
  }
  for (std::map<std::string, Long64_t>::const_iterator it =
         invalidByProcess.begin();
       it != invalidByProcess.end(); ++it) {
    std::cout << "  " << it->first << ": " << it->second
              << std::endl;
  }

  std::cout << "\ncreatorProcessName for valid rows" << std::endl;
  if (validByProcess.empty()) {
    std::cout << "  none" << std::endl;
  }
  for (std::map<std::string, Long64_t>::const_iterator it =
         validByProcess.begin();
       it != validByProcess.end(); ++it) {
    std::cout << "  " << it->first << ": " << it->second
              << std::endl;
  }

  std::cout << "\ncreatorProcessName for valid rows with "
            << "primaryGammaTrackID == -1" << std::endl;
  if (validPrimaryUnknownByProcess.empty()) {
    std::cout << "  none" << std::endl;
  }
  for (std::map<std::string, Long64_t>::const_iterator it =
         validPrimaryUnknownByProcess.begin();
       it != validPrimaryUnknownByProcess.end(); ++it) {
    std::cout << "  " << it->first << ": " << it->second
              << std::endl;
  }

  std::cout << "\nParent gamma trackIDs for valid rows with "
            << "primaryGammaTrackID == -1" << std::endl;
  if (validPrimaryUnknownByParentGamma.empty()) {
    std::cout << "  none" << std::endl;
  }
  for (std::map<int, Long64_t>::const_iterator it =
         validPrimaryUnknownByParentGamma.begin();
       it != validPrimaryUnknownByParentGamma.end(); ++it) {
    std::cout << "  parentGammaTrackID=" << it->first
              << ": " << it->second << std::endl;
  }

  TTree* gammaTree = nullptr;
  TTree* gammaEntryTree = nullptr;
  file->GetObject("GammaInteractionTree", gammaTree);
  file->GetObject("GammaMcpEntryTree", gammaEntryTree);

  std::set<GammaKey> gammaTrackKeys;
  std::map<GammaKey, std::set<int> > parentIDsByGamma;
  std::map<GammaKey, std::vector<double> > energiesByGamma;
  LoadGammaTrackInfo(gammaTree,
                     gammaTrackKeys,
                     parentIDsByGamma,
                     energiesByGamma);
  LoadGammaTrackInfo(gammaEntryTree,
                     gammaTrackKeys,
                     parentIDsByGamma,
                     energiesByGamma);

  if (!unknownPrimaryGammaKeys.empty()) {
    Long64_t uniqueGammaCount = 0;
    Long64_t gammaParentAppearsAsGamma = 0;
    Long64_t gammaParentDoesNotAppearAsGamma = 0;
    Long64_t gammaParentUnknownNoInfo = 0;
    Long64_t electronRowsParentAppearsAsGamma = 0;
    Long64_t electronRowsParentDoesNotAppearAsGamma = 0;
    Long64_t electronRowsParentUnknownNoInfo = 0;
    std::map<int, Long64_t> gammaParentIDDistribution;
    std::map<int, Long64_t> energyBins100keV;
    Double_t minEnergy = 1.0e99;
    Double_t maxEnergy = -1.0e99;
    Long64_t energySampleCount = 0;

    for (std::set<GammaKey>::const_iterator keyIt =
           unknownPrimaryGammaKeys.begin();
         keyIt != unknownPrimaryGammaKeys.end(); ++keyIt) {
      ++uniqueGammaCount;

      const Long64_t electronRowsForGamma =
        unknownPrimaryRowsByGamma[*keyIt];
      const std::map<GammaKey, std::set<int> >::const_iterator
        parentFound = parentIDsByGamma.find(*keyIt);

      if (parentFound == parentIDsByGamma.end() ||
          parentFound->second.empty()) {
        ++gammaParentUnknownNoInfo;
        electronRowsParentUnknownNoInfo += electronRowsForGamma;
      } else {
        bool anyParentAppearsAsGamma = false;
        for (std::set<int>::const_iterator parentIt =
               parentFound->second.begin();
             parentIt != parentFound->second.end(); ++parentIt) {
          ++gammaParentIDDistribution[*parentIt];
          if (gammaTrackKeys.find(GammaKey(keyIt->first, *parentIt)) !=
              gammaTrackKeys.end()) {
            anyParentAppearsAsGamma = true;
          }
        }

        if (anyParentAppearsAsGamma) {
          ++gammaParentAppearsAsGamma;
          electronRowsParentAppearsAsGamma += electronRowsForGamma;
        } else {
          ++gammaParentDoesNotAppearAsGamma;
          electronRowsParentDoesNotAppearAsGamma += electronRowsForGamma;
        }
      }

      const std::map<GammaKey, std::vector<double> >::const_iterator
        energiesFound = energiesByGamma.find(*keyIt);
      if (energiesFound != energiesByGamma.end()) {
        for (std::vector<double>::const_iterator energyIt =
               energiesFound->second.begin();
             energyIt != energiesFound->second.end(); ++energyIt) {
          const Double_t energy = *energyIt;
          if (energy < minEnergy) {
            minEnergy = energy;
          }
          if (energy > maxEnergy) {
            maxEnergy = energy;
          }
          ++energyBins100keV[static_cast<int>(energy/100.0)];
          ++energySampleCount;
        }
      }
    }

    std::cout << "\nParent-chain diagnostic for valid rows with "
              << "primaryGammaTrackID == -1" << std::endl;
    std::cout << "  electron rows with primaryGammaTrackID == -1: "
              << validPrimaryUnknown << std::endl;
    std::cout << "  unique incriminated gamma tracks: "
              << uniqueGammaCount << std::endl;
    std::cout << "  unique gammas whose parentID appears as gamma track: "
              << gammaParentAppearsAsGamma << std::endl;
    std::cout << "  unique gammas whose parentID does not appear as gamma: "
              << gammaParentDoesNotAppearAsGamma << std::endl;
    std::cout << "  unique gammas with no parent info in gamma trees: "
              << gammaParentUnknownNoInfo << std::endl;
    std::cout << "  electron rows in parent-appears category: "
              << electronRowsParentAppearsAsGamma << std::endl;
    std::cout << "  electron rows in parent-not-gamma category: "
              << electronRowsParentDoesNotAppearAsGamma << std::endl;
    std::cout << "  electron rows with no parent info: "
              << electronRowsParentUnknownNoInfo << std::endl;

    std::cout << "\nparentID distribution for incriminated gammas"
              << std::endl;
    if (gammaParentIDDistribution.empty()) {
      std::cout << "  none" << std::endl;
    }
    for (std::map<int, Long64_t>::const_iterator it =
           gammaParentIDDistribution.begin();
         it != gammaParentIDDistribution.end(); ++it) {
      std::cout << "  parentID=" << it->first
                << ": " << it->second << " unique gamma(s)"
                << std::endl;
    }

    std::cout << "\nenergy distribution for incriminated gammas"
              << std::endl;
    if (energySampleCount == 0) {
      std::cout << "  none" << std::endl;
    } else {
      std::cout << "  energy samples from gamma trees: "
                << energySampleCount << std::endl;
      std::cout << "  energy range [keV]: "
                << minEnergy << " -> " << maxEnergy << std::endl;
      for (std::map<int, Long64_t>::const_iterator it =
             energyBins100keV.begin();
           it != energyBins100keV.end(); ++it) {
        const int low = 100*it->first;
        const int high = low + 100;
        std::cout << "  [" << low << ", " << high
                  << ") keV: " << it->second << std::endl;
      }
    }
  }

  std::cout << "\nConvention reminder:" << std::endl;
  std::cout << "  primaryGammaTrackID is the primary generated gamma "
            << "when the gamma ancestry is known through gamma parents."
            << std::endl;
  std::cout << "  If a gamma is produced by a non-gamma parent, the "
            << "current safe convention is primaryGammaTrackID = -1."
            << std::endl;

  file->Close();
  delete file;
}
