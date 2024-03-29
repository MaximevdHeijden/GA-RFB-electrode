# Example script 
To run an example script, the "GA_main_Windows" script from this folder should be running without changing any parameters, resulting in numerous ".pnm" output files and one ".xlsx" file. The results of this run should be similar to the data stored in the "Genetic_Algorithm/GA_2/Reference" folder, where the ".pnm" files are the best individuals of the first and last generation determined from the "Cubic_Reference.xlsx" file.\
The pore size distribution, pore and throat size information, polarization curve, vtk files, and velocity profiles that can be obtained from the ".pnm" files can be found in the "output/GA" folder and can be obtained by running the "Network_properties_Windows" and "GA_to_VTK_and_properties_Windows" scripts. More information about the output of the scripts and the data visualization can be found in this [BLOG](https://maximevdheijden.github.io/2023/12/15/VisualizationParaview/).

## The following data could be obtained from the example script:
From the "GA_to_VTK_and_properties_Windows" script using the output "..._velocity.vtp" file, the pore sizes and for example the velocity, concentration, and pore current profiles could be visualized:\
<img src="output/Example_1.png" alt="color photo ftl" width="50%" height="auto" />\
*Figure 1: The geometrical evolution of the reference system. The networks of the first and final (1000th) generations are shown, displaying the: (a) pore diameter evolution, (b) the throat absolute velocity, (c) the pore concentration, and (d) the absolute pore current. With the flow in the y-direction and the thickness in the z-direction with the membrane facing to the top.*

Using the output ".xlsx" file obtained from "GA_main_Windows", the fitness values, electrical power, and pumping power can be plotted over the generations. The polarization curves can be obtained from the ".xlsx" output file of "GA_to_VTK_and_properties_Windows".\
<img src="output/Example_2.png" alt="color photo ftl" width="50%" height="auto" />\
*Figure 2: The results of the genetic algorithm optimized for the VO2+/VO2+ electrolyte, evaluated for the reference system (cubic network with mutation and a flow-through flow field), with: (a) the structure evolution over 1000 generations with the flow in the y-direction and the thickness in the z-direction with the membrane facing to the top, (b) the maximum fitness evolution, (c) the maximum electrical power and pumping power evolution, and (d) the polarization curve for the first and last generation.*

Finally, the pore and throat size distributions can be obtained from the ".xlsx" file generated using the "Network_properties_Windows" script.\
<img src="output/Example_3.png" alt="color photo ftl" width="50%" height="auto" />\
*Figure 3: The (a) pore and (b) throat size distributions of the first and final generation, showing the pore or throat count and cumulative normalized volume distributions (divided in 2 μm pore or throat sized bins) for the reference system.*
