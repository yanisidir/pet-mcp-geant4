# Minimal Geant4 Example

This folder is a clean, self-contained Geant4 starting point. It does not use
ROOT, the old MCP geometry, parameterisations, custom sensitive detectors, or
custom physics.

## Project Structure

```text
minimal_example/
  CMakeLists.txt
  minimal_geant4.cc
  include/
    ActionInitialization.hh
    DetectorConstruction.hh
    EventAction.hh
    PrimaryGeneratorAction.hh
    RunAction.hh
    SteppingAction.hh
  src/
    ActionInitialization.cc
    DetectorConstruction.cc
    EventAction.cc
    PrimaryGeneratorAction.cc
    RunAction.cc
    SteppingAction.cc
  macros/
    run.mac
    vis.mac
  README_MINIMAL.md
```

## Role of Each Geant4 Class

- `minimal_geant4.cc`: program entry point. It creates the run manager, geometry,
  physics list, user actions, and optional visualization. It uses one thread so
  the console output stays easy to read.
- `DetectorConstruction`: defines the geometry and materials. The geometry is one
  air world box and one silicon detector box.
- `PrimaryGeneratorAction`: creates one primary particle per event using a
  `G4ParticleGun`.
- `ActionInitialization`: registers the primary generator, run action, event
  action, and stepping action.
- `RunAction`: prints when a run starts and ends.
- `EventAction`: prints when each event starts and ends.
- `SteppingAction`: prints when a particle crosses from the world into the
  detector volume.
- `FTFP_BERT`: the standard Geant4 reference physics list used by the main
  program.

## Compile

From this repository:

```bash
cd minimal_example
cmake -S . -B build
cmake --build build
```

If your Geant4 environment is not already loaded, source your Geant4 setup file
first. For example:

```bash
source /path/to/geant4-install/bin/geant4.sh
```

## Run in Batch Mode

```bash
cd minimal_example/build
./minimal_geant4 macros/run.mac
```

The batch macro runs five events with a 1 MeV gamma starting at
`(0, 0, -40 cm)` and moving along `+z`.

## Run With Visualization

```bash
cd minimal_example/build
./minimal_geant4
```

This opens an interactive Geant4 UI session and executes `macros/vis.mac`.

## Change the Primary Particle

Edit `macros/run.mac` or type these commands in the interactive UI after
`/run/initialize` and before `/run/beamOn`:

```text
/minimal/gun/particle e-
/minimal/gun/energy 5 MeV
/minimal/gun/position 0 0 -40 cm
/minimal/gun/direction 0 0 1
```

Useful particle names include `gamma`, `e-`, `e+`, `proton`, and `mu-`.

The direction vector does not need to be normalized in the macro. The code calls
`.unit()` before passing it to the particle gun.

## What Happens During One Event

1. `PrimaryGeneratorAction` creates one primary particle.
2. Geant4 transports the particle through the air world using `FTFP_BERT`.
3. If the particle crosses into the silicon detector box, `SteppingAction`
   prints the particle name, track ID, kinetic energy, and entry position.
4. `EventAction` prints that the event has ended.

## Add ROOT Output Later

A simple way to add ROOT back is:

1. Add `find_package(ROOT REQUIRED COMPONENTS RIO Tree)` in `CMakeLists.txt`.
2. Create a small output class that owns a `TFile` and `TTree`.
3. Fill the tree from `SteppingAction` or `EventAction`.
4. Open the file in `RunAction::BeginOfRunAction`.
5. Write and close it in `RunAction::EndOfRunAction`.

Keep ROOT isolated in one or two files so the geometry and generator remain easy
to understand.

## Add MCP Geometry Later

Rebuild the MCP geometry inside `DetectorConstruction` in small steps:

1. Keep the current world volume.
2. Add one MCP component at a time as a named logical volume.
3. Verify each new component in visualization before adding the next one.
4. Add parameterisation only after a manually placed version is understood.
5. Add sensitive detector logic only after the passive geometry is stable.

This keeps the minimal example useful as a known-good base while advanced
features are reintroduced.
