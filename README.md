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

## Build

```bash
mkdir -p build
cd build
cmake ..
make -j
```

## Run

```bash
cd build
./minimal_geant4 macros/run.mac
```

When the simulation is launched from `build/`, the ROOT output is written as:

```text
build/mcp_output.root
```

In the current serial configuration, the normal output file is
`mcp_output.root`. In multithreaded configurations, if enabled, ROOT files may
be written per worker thread, for example `mcp_output_t0.root`,
`mcp_output_t1.root`, etc.

## Visualization

```bash
cd build
./minimal_geant4 macros/vis_current.mac
```

Some `vis_*` macros are kept for historical views of older multi-MCP
geometries. If a macro references volumes such as `MCP_plus_body_1` while the
current code uses `kNumberOfMCPs = 1`, prefer `vis_current.mac`.

```bash
./minimal_geant4 macros/vis_global.mac
./minimal_geant4 macros/vis_plus_z.mac
./minimal_geant4 macros/vis_minus_z.mac
./minimal_geant4 macros/vis_face_plus_z.mac
./minimal_geant4 macros/vis_export.mac
```

## ROOT Analysis

From the project directory:

```bash
cd /path/to/pet-mcp-geant4
root -l -b -q 'script/inspect_electrons.C("build/mcp_output.root")'
```

Most analysis macros default to reading `build/mcp_output.root`, so launch them
from the project root unless you pass an explicit ROOT filename.

Figures produced by ROOT analysis macros are written to:

```text
Fig/current/
```

Older figures are archived in `Fig/old_*` directories.

Generated build files, ROOT outputs, and visualization exports are excluded
from Git.
