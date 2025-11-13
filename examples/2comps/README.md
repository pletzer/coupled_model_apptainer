# Two-component coupled example
i

This example shows how to create a 2-component coupled model. 
Each component lives on a different grid and data are 
automatically regridded when imported. 

|----------------
| component | import  | export |
| :-------- | :------ | :------|
| atm       | sst     | pmsl, rsns |
| ocn       | pmsl, rsns | sst  |


## How to build the code

```
mkdir build
cd build
apptainer shell /nesi/nobackup/nesi99999/pletzera/sifs/esmfenv.sif
cmake ..
make
```

## How to run the code

```
mpiexec -n 4 ./2com_time_example
```

