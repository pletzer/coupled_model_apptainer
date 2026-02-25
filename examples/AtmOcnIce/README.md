# Concurrent atm–ocn–ice Coupling Example

This example demonstrates how to couple three model components:

- **atm** (atmosphere)  
- **ocn** (ocean)  
- **ice** (sea ice)

The components run **concurrently**, each on its own set of PETs.

---

## Model Architecture

Each component:

- Runs on a dedicated subset of PETs  
- Owns its own 2D structured grid  
- Uses its own regular MPI domain decomposition  

Field and flux exchanges are handled through a **connector component**, which performs conservative regridding between the different component grids.

Currently, a single connector is used because all exchanges use conservative remapping.

---

## Building the Code

```bash
mkdir build
cd build
apptainer --bind /nesi/nobackup/ shell /nesi/nobackup/nesi99999/pletzera/sifs/esmfenv86.sif
cmake ..
make
```
This assumes the appropriate ESMF-enabled container environment.

## Running the Code

```
mpiexec -n 6 ./esmApp
```
The domain decompositions are currently hardwired.
You must run on exactly 6 MPI processes.

## Code Structure and Extension Guide

### Driver: esmfApp.F90
This is the main driver.
Responsibilities:

 * Defines the coupling configuration
 * Registers exchanged fields
 * Manages the NUOPC field dictionary

All exchanged fields must have valid standard names.
The default standard name list is limited, so new entries must be added manually.
Example:

```
call NUOPC_FieldDictionaryAddEntry( &
    standardname="ice_albedo", &
    canonicalunits="1", &
    rc=rc)
```

### Component Setup: esmf.F90

This file instantiates all components and connectors
Defines PET allocations for each component
Controls component concurrency layout
To modify processor allocation, look for arrays such as:
`atmPetList`.

To add new connectors (e.g., ICE → OCN), edit this file.

### Connector: conn.F90

Implements field remapping between components.
Uses conservative regridding
Handles all inter-component exchanges
At present, a single connector supports all exchanges because they share the same remapping method.

### Model Components

 * atm.F90
 * ocn.F90
 * ice.F90

Each component:
 * Creates its own grid
 * Defines import/export fields
 * Initializes and advances exported fields
 * Receives imported fields via NUOPC
 * Fields can be written to VTK format for visualization:
`call write_vtk(field_sst, filename)`

## Notes and Limitations
 * Fixed MPI layout (6 ranks required)
 * Structured grids only
 * Conservative interpolation only
 * Minimal field dictionary (manual extensions required)








