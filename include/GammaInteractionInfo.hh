#ifndef GAMMA_INTERACTION_INFO_HH
#define GAMMA_INTERACTION_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Structure décrivant une interaction gamma importante
// observée pendant le stepping.
//
// Une instance = une interaction gamma enregistrée.
//
// Exemples d'interactions sauvegardées :
// - effet photoélectrique (phot)
// - diffusion Compton (compt)
// - conversion de paire (conv)
// - diffusion Rayleigh (Rayl)
//
// Cette structure est utilisée pour :
// - stocker les interactions dans EventAction ;
// - remplir le TTree "gamma_interactions" ;
// - analyser le comportement des photons dans le MCP.
struct GammaInteractionInfo
{
  // Numéro de l'événement Geant4.
  G4int eventID;

  // Identifiant de la track gamma concernée.
  G4int trackID;

  // Identifiant de la particule mère.
  //
  // Pour un gamma primaire :
  // parentID = 0
  G4int parentID;

  // Identifiant numérique du processus.
  //
  // 1 = phot
  // 2 = compt
  // 3 = conv
  // 4 = Rayl
  G4int processID;

  // Flag indiquant si cette interaction correspond
  // à la première interaction du gamma primaire.
  //
  // 0 = non
  // 1 = oui
  G4int isFirstPrimaryInteraction;

  // Nom du processus responsable de l'interaction.
  //
  // Exemples :
  // "phot"
  // "compt"
  // "conv"
  // "Rayl"
  G4String processName;

  // Nom du volume dans lequel l'interaction a eu lieu.
  //
  // Exemples :
  // "MCP_body"
  // "MCP_channel"
  G4String volumeName;

  // Energie déposée pendant ce step.
  //
  // Pour certaines interactions gamma,
  // cette valeur peut être faible voire nulle.
  G4double energyDeposit;

  // Energie cinétique du gamma juste avant l'interaction.
  G4double kineticEnergyBefore;

  // Energie cinétique du gamma juste après l'interaction.
  G4double kineticEnergyAfter;

  // Longueur du step associé à cette interaction.
  G4double stepLength;

  // Temps global auquel l'interaction s'est produite.
  G4double globalTime;

  // Position de l'interaction.
  //
  // Dans cette simulation, elle correspond généralement
  // au PostStepPoint où le processus a été déclenché.
  G4ThreeVector position;

  // Direction du gamma après l'interaction.
  //
  // Particulièrement utile pour :
  // - les diffusions Compton ;
  // - les diffusions Rayleigh ;
  // - l'étude de la trajectoire du photon dans le MCP.
  G4ThreeVector directionAfter;
};

#endif