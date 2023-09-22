### Input Files 

This folder has input files required to run simulations for all water models described in the paper. 

* `Lsolv_template.top` is a general template topology file. The values
  required to complete the template can be found in
  `Force_field_parameters.csv` file.

* A sample initial .gro file with 521 atoms and box size 2.6 nm om
  each side is provided, as well as min.mdp, eq.mdp and run.mdp file
  which were used to run the GROMACS md simulation.

* `write_tops_from_template.py` is a python script that takes `Force_field_parameters.csv` and populates
  the individual force field files from the `Lsolv_template.top`