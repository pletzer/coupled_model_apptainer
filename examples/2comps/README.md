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
apptainer shell /nesi/nobackup/nesi99999/pletzera/sifs/esmfenv86.sif
cmake ..
make
```

## How to run the code

Under Apptainer:
```
mpiexec -n 4 ./2comp_time_example
```
This will generate `PET*.ESMF_LogFile` files.
