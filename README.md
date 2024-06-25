# Monte Carlo sampling of fluxes for finite temperature calculation of the Kitaev model


## Usage: 
To compile and set the correct enviroment variables, in the main directory of the code,
```bash
 $  source build.sh
```
To run the code, you have to be in a diferent directory with a parameter file like the one in "scripts/Start". Then to run the code, write
```bash
 $  $KITAEV_DIR/src/main
```
To analyse the results, 
```bash
 $  $KITAEV_DIR/src/ana obs_1 obs_2 ...
```
To analyse all of the observables, 
```bash
 $  $KITAEV_DIR/src/ana *
```

## Requirements:
The following python packages are required to run the program:
* C compiler
* F90 compiler
* GSL (GNU Scientific Library)





