# Two-component coupled example


This example shows how to create a 2-component coupled model. 
Each component lives on a different grid and data are 
automatically regridded onto the target component's grid. 

The data exchanged are sea surface temperature, air pressure at sea_level (pmsl) and
surface net downward shortwave flux (rsns). 

The coupled system advances from 2010-06-01:00:00 to 2010-06-01:01:00 with coupling occurring every 15 minutes. 


| component | import     | export     |
| :-------- | :------    | :------    |
| atm       | sst        | pmsl, rsns |
| ocn       | pmsl, rsns | sst        |


In addition to atm and ocn, there is a third component (esm), which handles the coupling.


## How to build the code

```
mkdir build
cd build
apptainer shell -B /nesi/nobackup /nesi/nobackup/nesi99999/pletzera/sifs/esmfenv86.sif
cmake ..
make
```

## How to run the code

Under Apptainer:
```
mpiexec -n 4 ./2comp_time_example
```
This will generate `PET*.ESMF_LogFile` files.

---

## Files Overview

### **1. `esmfApp.F90`**
- Main program that initializes and runs the ESM application.
- Responsibilities:
  - Call `ESMF_Initialize` and `ESMF_Finalize`.
  - Create the top-level ESM component (`esmComp`) using `ESMF_GridCompCreate`.
  - Set the NUOPC `SetServices` routine for the component.
  - Call `ESMF_GridCompInitialize`, `ESMF_GridCompRun`, and `ESMF_GridCompFinalize`.
  - Destroy the ESM component.
- Acts as the **entry point** of the application.

---

### **2. `esmf.F90`**
- Module that specializes the generic ESM component.
- Responsibilities:
  - Provides the `SetServices` routine for the driver component.
  - Derives the driver from `NUOPC_Driver`.
  - Specializes the driver with the `SetModelServices` routine.
  - Adds child components (ATM and OCN) to the driver.
  - Adds connectors between ATM → OCN and OCN → ATM.
  - Sets up the application clock with start/stop time and time step.

---

### **3. `atm.F90`**
- Atmospheric model component.
- Responsibilities:
  - Provides `SetServices` routine to derive from `NUOPC_Model` and specialize `Advertise`, `RealizeProvided`, and `Advance` routines.
  - **Advertise**:
    - Importable field: `sea_surface_temperature` (sst)
    - Exportable fields: `air_pressure_at_sea_level` (pmsl), `surface_net_downward_shortwave_flux` (rsns)
  - **Realize**:
    - Reads configuration from `2comp_time_example.nml` for grid dimensions and coordinates.
    - Creates `ESMF_Grid` objects for input and output fields.
    - Creates and realizes the importable and exportable fields.
  - **Advance**:
    - Moves the model forward in time according to its internal clock.
    - Logs current and next time steps.

---

### **4. `ocn.F90`**
- Ocean model component.
- Responsibilities:
  - Provides `SetServices` routine to derive from `NUOPC_Model` and specialize `Advertise`, `RealizeProvided`, `SetClock`, and `Advance` routines.
  - **Advertise**:
    - Importable fields: `air_pressure_at_sea_level` (pmsl), `surface_net_downward_shortwave_flux` (rsns)
    - Exportable field: `sea_surface_temperature` (sst)
  - **Realize**:
    - Reads configuration from `2comp_time_example.nml`.
    - Creates `ESMF_Grid` objects for input and output fields.
    - Creates and realizes importable and exportable fields.
  - **SetClock**:
    - Initializes the internal clock with a stability time step (5 minutes in the example).
  - **Advance**:
    - Moves the ocean model forward in time according to its internal clock.
    - Logs current and next time steps.

---

## Summary of Data Flow

- ATM exports `air_pressure_at_sea_level` and `surface_net_downward_shortwave_flux`.
- OCN imports these fields from ATM.
- OCN exports `sea_surface_temperature`.
- ATM imports `sea_surface_temperature` from OCN.
- Connector components (`cplSS`) handle the mapping of fields between ATM and OCN grids.

---

## Accessing Imported Fields

To access an imported field within a component (e.g., ATM):

```fortran
! Local variables
type(ESMF_State) :: importState
type(ESMF_Field) :: sstField
integer :: rc

! Query the import state
call NUOPC_ModelGet(model, importState=importState, rc=rc)

! Retrieve the imported field by name
call ESMF_StateGet(importState, field=sstField, itemName="sst", rc=rc)

! Access the data array
call ESMF_FieldGet(sstField, farrayPtr=dataPtr, rc=rc)
! dataPtr now points to the Fortran array of field values
```

itemName corresponds to the name given in Advertise.
The retrieved pointer (dataPtr) can then be used for calculations or updates in your component.
The connector automatically handles interpolation/remapping between source and local grids.

---

## Field Remapping (Regridding)

Remapping occurs when fields are transferred between components with different grids. In this example, the connectors handle this automatically.

### Remap Methods

 * CONSERVE: Conservative remapping. Preserves integrals (e.g., total mass/energy). Best for fluxes.
 * BILINEAR: Bilinear interpolation. Smooth but may not conserve fluxes.
 * PATCH: Higher-order patch interpolation.
 * NEAREST_STOD: Nearest-neighbor interpolation from source to destination.

### Setting the Remap Method in a Connector

```fortran
type(ESMF_CplComp) :: connector
integer :: rc

! Add connector ATM->OCN
call NUOPC_DriverAddComp(driver, srcCompLabel="ATM", dstCompLabel="OCN", &
    compSetServicesRoutine=cplSS, comp=connector, rc=rc)

! Set the remap method
call NUOPC_CompAttributeSet(connector, name="RemapMethod", value="CONSERVE", rc=rc)
```
You can switch between "BILINEAR" or "CONSERVE" depending on whether you want smooth interpolation or conservative fluxes.

The same connector can handle multiple fields. Each field will be regridded according to the same RemapMethod.
If you need different interpolation methods for different fields, you would either:
 * Create separate connectors per method, or
 * Use advanced regridding inside the connector’s Regrid object, which can handle per-field weights.


### How it Works at Runtime
 1. During Realize, the source and destination fields are created.
 2. The connector internally builds an ESMF_Regrid object using:
   * Source grid
   * Destination grid
   * Remap method
 3. At each time step, ESMF_Regrid is applied automatically, and the import field is updated.

### Advanced Options

 * Masking: skip land or ocean points.
 * Remap weights file: precompute weights for static grids for performance.
 * Higher-order conservative: RemapOrder=2 or higher.


## Notes

- All components are derived from NUOPC base classes (`NUOPC_Model`, `NUOPC_Driver`).
- Fields are created with `ESMF_Grid` objects specifying the spatial domain.
- The time step for each component can be controlled via the internal clock.
- You can enable/disable import fields or connectors by commenting/uncommenting the corresponding macros (`WITHIMPORTFIELDS`, `WITHCONNECTORS`).

