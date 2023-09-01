# Stars Collisions  :collision:

[![Up to Date](https://github.com/ikatyang/emoji-cheat-sheet/workflows/Up%20to%20Date/badge.svg)](https://github.com/JuanManuelPacheco/Star_Collisions/actions)

This repository was created to record and save the progress on the final project of the computational astrophysics subject part of the Master in Astrophysics and Cosmology of the University of Padova during the 2021-2022 year.

The central goal of the project was to analyze a set of star collision simulations generated using the sph **StarSmasher** code and to determine:

1) From which instant in the time evolution of a collision the final remnant is considered to be sufficiently relaxed.

2) To analyze the rotational behavior of the remnant.

## Structure of this repository:

* **Bash:** This folder contains all the bash scripts generated in the development of this project. *script_ascii.sh* produce the ascii files for all the snapshots using **splash**, *script_Menc.sh* applies the *Menc_Munb.py* module on all the ascii files previously produced, *script_L.sh* applies the *L_evo.py* module and finally  *script_Main.sh* applies the *Main.py* module.  

* **Data:** This folder contains the main .txt files produced by the aplication of the bash scripts over all snapshots divided by folders of each simulation and the video for the evolution during the collisions. The file *legend.txt* contains the explanation of all the main characteristics of each collision simulated.

* **References:** This folder contains all the pdf's of the references used in the development of this project.

* **Python:** This folder contains all the .py files generated in the development of this project. *Menc_Munb.py* computes the mass profiles for all the snapshots, *L_evo.py* calculates the radial, tangential and vertical velocity of each simulation, and *Main.py* module performs all these calculations on one shot.  

* **Stars_collisions:** This jupyter notebook was generated as the report needed in the presentation of this project. It contains all the explanations and documentation on the functions and procedures needed to achieve the main objectives.
