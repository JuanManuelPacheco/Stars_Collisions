# Stars Collisions  :collision:

[![Up to Date](https://github.com/ikatyang/emoji-cheat-sheet/workflows/Up%20to%20Date/badge.svg)](https://github.com/JuanManuelPacheco/Star_Collisions/actions)

This repository was created to record and save the progress on the thesis project part of the Master in Astrophysics and Cosmology of the University of Padova during the 2021-2023 year.

* The goal of this master thesis is to follow the hybrid methodology, proposed by *Ballone, et al. (2023)* and *Costa, et al. (2022)*, to explore other types of stellar collisions by the analisis of a set of star collision simulations generated using the SPH **StarSmasher** code. The idea is to vary the evolutionary stages of both stars before the encounter and the orbital parameters of the collision, identifying the influence they have on the post-coalescing stellar structure. Finally, the parameter space of the stellar collisions formation channel is expected to be explored in more detail throughout this thesis thanks to the generation and analysis of a set of stellar profiles and hydrodynamical collisions, concluding with the necessary characteristics to populate the pair-instability mass gap with black holes generated through the encounter of massive stars in young stellar clusters.

## Structure of this repository:

* **Bash:** This folder contains all the bash scripts generated in the development of this project. *script_ascii.sh* produce the ascii files for all the snapshots using **splash**, *script_assign.sh* applies the *Assignation.py* module, *script_E.sh* applies the *Energy_Computation.py* module,  *script_entropy.sh* applies the *Entropy.py* module, *script_Menc.sh* applies the *Menc_Munb.py* module on all the ascii files previously produced, *script_L.sh* applies the *L_evo.py* module and finally *script_Main.sh* applies the *Main.py* module.  

* **Data:** This folder contains the main .txt files produced by the aplication of the bash scripts over all snapshots divided by folders of each simulation. The file *coll_legend.txt* contains the explanation of all the main characteristics of each collision simulated. The MESA and PARSEC folders contains the estellar evolution profiles of the relaxed stars.

* **References:** This folder contains all the pdf's of the references used in the development of this thesis.

* **Figures:** This folder contains all the png files included in the thesis document.

* **Python:** This folder contains all the .py files generated in the development of this project. *Assignation.py* associates each particle of the final collisional snapshot to its respective initial star and calculates the fraction of hydrogen and helium corresponding to it, *Energy_Computation.py* calculates the total and approximate energy of each particle sph for comparison, *Entropy.py* calculates the entropic variable of each of the particles in the desired snapshot,  *Menc_Munb.py* computes the mass profiles for all the snapshots, *L_evo.py* calculates the radial, tangential and vertical velocity of each simulation, and *Main.py* module performs all these calculations on one shot.  

* **Notebooks:** This folder contains all the  jupyter notebook generated for this project.
  
* **Master_Thesis.pdf:** This is the final pdf file submitted for the evalutation of the committee to obatin the Master's degree.
