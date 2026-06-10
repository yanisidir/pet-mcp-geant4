#ifndef MCP_PLATE_STATS_INFO_HH
#define MCP_PLATE_STATS_INFO_HH

#include "globals.hh"

// Une ligne par événement, côté et plaque MCP.
struct McpPlateStatsInfo
{
  G4int eventID;
  G4int side;
  G4int mcpIndex;
  G4int electronChannelCount;
};

#endif
