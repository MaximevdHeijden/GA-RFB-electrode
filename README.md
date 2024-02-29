# Genetic Algorithm scripts
## About
This repository contains a novel approach for the optimization of porous electrode microstructures from the bottom-up is presented that couples a genetic algorithm with a validated [pore network model](https://github.com/MaximevdHeijden/PNM-RFB-electrode).\
In our developed genetic algorithm (GA), optimized porous electrodes are computed from the bottom-up by natural selection of the fittest structures driven by the [theory of evolution](https://maximevdheijden.github.io/2023/11/23/Computationaloptimizationstrategies/). A GA has a high degree of evolutionary freedom and can easily be adapted and expanded depending on the operating conditions (potential, flow rate), network properties (artificial, extracted), system configuration (flow field, electrolyte), and desired optimization functions (mutation, migration). The framework was built to speed up and guide the design of electrodes for electrochemical applications, with a focus on redox flow batteries.\
GA-RFB-electrode is an open-source code written in Python using the open-source software OpenPNM.

## Installation
The scripts and functions in this folder use OpenPNM version 2.6.0. which can be installed using [OpenPNM documentation](https://openpnm.org/installation.html). Before running the code, minor adjustments need to be made to the OpenPNM documentation, which can be found in the “READ Me – OpenPNM changes” file. For the installation of OpenPNM from scratch follow the steps below:

1.	Download [Anaconda](https://www.anaconda.com/download/) and follow the installment procedure. The open-source Anaconda Distribution is the easiest way to perform Python/R data science and machine learning on Linux, Windows, and Mac OS X. Python is used as the programming language in this project and Spyder can be used as editor.

2.    Check pip installation (pip is the package manager of python) and install pipenv. First, open 'Anaconda prompt' via the Windows start menu. Check if pip is installed using:

      ```
      python -m pip --version
      ```

      If pip is installed, install pipenv using:

      ```
      python -m pip install --user pipenv
      ```

      You can get a warning that the installed packages are not on your PATH. Copy this path, go to your windows start menu and type in “Edit the system environment variables”, select “environment variables”, click on Path and select “Edit”. Then click “new” and paste the path you just copied.

3.    Install the latest version of [OpenPNM](https://openpnm.org/installation.html), following the “hard (but correct way)” will allow you to install the developing code:\
      •	First make a repository (folder) in which you want to save the code and copy the directory address.\
      •	Open Anaconda Prompt and go to the folder by using the "cd" command.\
      •	Clone the repo via the command:
               
            git clone https://github.com/PMEAL/OpenPNM
                  
      •	Enter the root folder of openPNM by using the "cd OpenPNM" command.\
      •	Enter the following command and note the space between the dot and e is not a typo:
            
            pip install --no-deps -e .
            
      •	Enter the following commands to install the openpnm dependencies:
       
            conda install --file requirements/conda.txt -c conda-forge
            pip install -e .

      •	Switch to version 2.6.0 of OpenPNM using:

            git checkout v2.6.0

5.    You need to install some dependencies using for example pip by entering the following command in AnacondaPrompt:
 
            pip install -r requirements.txt

6.    Now the codes should work.

## Documentation
This repository contains several scripts that will be used when running the code or for post-processing, including:

**GA_main:**
This script is used for running the algorithm and in which the desired operating conditions, network properties etc. must be selected.
Two versions exist of this code:
1.	**GA_main_Linux:** This script is suitable for running at Linux and contains a parallelized part, for which the number of cores should be defined. In the current script, the number of cores used is equal to the population size. 
2.	**GA_main_Windows:** This script is suitable for running at Windows and does not contain a parallelized part.
The algorithm only optimizes for the network size selected and not for the entire electrode length and only for one applied potential.

**GA_functions:**
The functions used in the algorithm are defined in this script, excluding the functions used in the integrated pore network model.
Two versions exist of this code:
1.	**GA_main_Linux:** This script is suitable for running at Linux.
2.	**GA_main_Windows:** This script is suitable for running at Windows.

**CustomFunctionsGA:**
The functions used in the pore network model, which is integrated in the genetic algorithm, are defined in this script. 
Two versions exist of this code:
1.	**CustomFunctionsGA:** This script is suitable for artificially generated networks such as the Cubic and Voronoi networks.
2.	**CustomFunctionsGA_FH23:** This script is suitable for tomographic extracted networks (the xyz coordinates in the provided extracted networks are different than for the cubic networks, affecting the boundary conditions).

**inputDict:**
The input parameters used in the algorithm are defined in this script. 
Four versions exist of this code:
1.	**inputDict_GA_V:** This script is suitable for artificially generated networks such as the Cubic and Voronoi networks and runs the algorithm for the vanadium chemistry.
2.	**inputDict_GA_V_FH23:** This script is suitable for tomographic extracted networks and runs the algorithm for the vanadium chemistry (the xyz coordinates in the provided extracted networks are different than for the cubic networks, affecting the boundary conditions).
3.	**inputDictTEMPO:** This script is suitable for artificially generated networks such as the Cubic and Voronoi networks and runs the algorithm for the TEMPO chemistry.
4.	**inputDictTEMPO_FH23:** This script is suitable for tomographic extracted networks and runs the algorithm for the TEMPO chemistry (the xyz coordinates in the provided extracted networks are different than for the cubic networks, affecting the boundary conditions).

**GA_to_VTK_and_properties_Windows:**
This script converts the networks generated and saved in ‘GA_main’ to VTK files which can be visualized in [Paraview](https://www.paraview.org/) (more information can be found in the following [BLOG](https://maximevdheijden.github.io/2023/12/15/VisualizationParaview/)). Moreover, the polarization plot will be obtained for a specified number of applied potentials, which can also be obtained using a network-in-series approach to mimic a larger electrode size.

**Network_properties_Windows:**
This script computes the pore size distribution of a network generated and saved in ‘GA_main’.

## Getting started
After installing OpenPNM version 2.6.0. and making the required adjustments, the genetic algorithm can be run for the desired operating conditions and network properties, which need to be specified in ‘GA_main’.\
The following properties can be adjusted in ‘GA_main’:\
•	Population size\
•	Number of parents\
•	Number of generations\
•	Mutation probability and range\
•	Initial guesses for the overpotential (should be linked to the applied overpotential, defined in the inputDict)\
•	Network properties including: the network shape, spacing, connectivity, and reference pore and throat diameters for cubic networks, the network shape and number of points for a voronoi network, and the network file for the extracted network.\
•	Merging and splitting probability and ratio.\
•	The network type: 0 = cubic network, 1 = extracted network, 2 = voronoi network.\
•	Merging and splitting option: 0 = only mutation, 1 = mutation and merging and splitting.\
•	Flow field type: 0 = flow-through flow field, 1 = interdigitated flow field.\
•	Reload data: option to continue the optimization at another time.\
•	Chemistry type: this repository contains the chemical information for the vanadium and TEMPO electrolytes.\
•	File name and folder.

More information regarding these parameters can be found in the publications mentioned under the ‘Cite’ section. **Care must be taken when changing the electrolyte and network types, as different inputDicts and customFunctions must be used. Moreover, the inputDict file used must additionally be changed in the customFunction script! Furthermore, the extracted network must be run at Linux as they don’t work yet on Windows.**\
**If other extracted networks are used, make sure you double check the xyz coordinates, as they affect the boundary conditions and thus the optimization.**\
**In addition, when using another electrolyte than vanadium, also adjust the end of the algorithm function in the ‘GA-functions script’, as the P_theory must include the full-cell potential of the desired electrolyte (which is 1.26 for the vanadium chemistry).**

Before running the code, the following folders need to be created:\
•	**output:** for saving the data from the 'GA_to_VTK_and_properties_Windows' script.\
•	**Genetic_Algorithm:** with a subfolder called **GA-2** for storing the data generated in the 'GA_main' script.\
•	The **input** folder in this repository additionally needs to be created when running the algorithm for extracted networks. Two examples, paper electrodes for a flow-through and an interdigitated flow field, are provided.

## Future alterations to the code
This first version of the code is not optimal in terms of user-friendliness and should therefore be polished further in the future including:\
•	The integration of the different chemistry and network versions.\
•	The saving of the networks: now all the networks in generations 1-10 and every 100th generation are saved instead of the best network in each generation.\
•	The selection of the mating pool will be changed as now only new individuals (after crossover) are created. In the newer version, the best-performing network in the previous generation (without crossover) will be used in the next generation as well. 

## Future research directions
The presented GA is a proof-of-concept and should be further developed in future work. Possible extensions or improvements of the modeling framework are provided in the following [BLOG](https://maximevdheijden.github.io/2023/11/22/GAprospects/).

## Contribution
GA-RFB-electrode is developed using a wide range of open-source tools, including OpenPNM. The code has initially been developed as part of a PhD research project and two master projects, but further contributions are very welcome.  

## Publications and citations
The code has been developed and described in the following two publications and PhD thesis. Please cite them when referring to the algorithm: 

```bash
@phdthesis{heijden2023engineering,
      title = {Engineering porous electrodes for redox flow batteries - modeling, diagnostics, and manufacturing},
      author = {Maxime van der Heijden},
      isbn = {978-90-386-5878-0},
      year = {2023},
      publisher = {Eindhoven University of Technology},
}
@article{heijden2023versatile,
      title = {A versatile optimization framework for porous electrode design},
      author = {Maxime van der Heijden and Gabor Szendrei and Victor de Haas and Antoni Forner-Cuenca},
      journal = {ChemRxiv},
      year = {2023},
      doi = {10.26434/chemrxiv-2023-nzr36},
}
@article{gorp2023bottomup,
      title = {Bottom-up design of porous electrodes by combining a genetic algorithm and a pore network model},
      author = {Rik van Gorp and Maxime van der Heijden and Mohammad Amin Sadeghi and Jeffrey Gostick and Antoni Forner-Cuenca},
      journal = {Chemical Engineering Journal},
      volume = {455},
      pages = {139947},
      year = {2023},
      issn = {1385-8947},
      doi = {10.1016/j.cej.2022.139947},
}
```

For referring to the [pore network model](https://github.com/MaximevdHeijden/PNM-RFB-electrode), please cite the following publications: 

```bash
@article{heijden2022assessing,
      title = {Assessing the versatility and robustness of pore network modeling to simulate redox flow battery performance},
      author = {Maxime van der Heijden and Rik van Gorp and Mohammad Amin Sadeghi and Jeffrey Gostick and Antoni Forner-Cuenca},
      journal = {Journal of the Electrochemical Society},
      volume = {169},
      pages = {040505},
      year = {2022},
      doi = {10.1149/1945-7111/ac5e46},
}
@article{munoz2023versatile,
      title = {Understanding the role of electrode thickness on redox flow cell performance},
      author = {Vanesa Muñoz-Perales and Maxime van der Heijden and Victor de Haas and Jacky Olinga and Marcos Vera and Antoni Forner-Cuenca},
      journal = {ChemElectroChem},
      year = {2023},
      doi = {10.1002/celc.202300380},
}
```
