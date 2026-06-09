# PET MCP Geant4

Geant4 simulation of a simplified PET detector using two symmetric MCP stacks.
Each event generates two back-to-back 511 keV gamma photons. The simulation
records electrons reaching MCP channels, photons leaving the detector, and
event-level coincidence statistics.

## Main Features

- MCP stacks on the `+z` and `-z` sides
- Parameterised MCP channels with configurable chevron angles
- Livermore electromagnetic physics
- Back-to-back 511 keV PET source
- Per-thread ROOT output for multithreaded runs
- Electron hit counts by side and MCP plate
- PET coincidence efficiency
- Geant4 visualization and export macros

## Requirements

- Geant4 11 with UI and visualization support
- ROOT
- CMake 3.16 or newer
- A C++11-compatible compiler

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Run

```bash
cd build
./minimal_geant4 macros/run.mac
```

In multithreaded mode, ROOT files are written per worker thread:

```text
mcp_output_t0.root
mcp_output_t1.root
...
```

## Visualization

```bash
cd build
./minimal_geant4 macros/vis_global.mac
./minimal_geant4 macros/vis_plus_z.mac
./minimal_geant4 macros/vis_minus_z.mac
./minimal_geant4 macros/vis_face_plus_z.mac
./minimal_geant4 macros/vis_export.mac
```

## ROOT Analysis

From the project directory:

```bash
root -l -b -q 'inspect_electrons.C("build/mcp_output_t0.root")'
```

Generated build files, ROOT outputs, and visualization exports are excluded
from Git.
