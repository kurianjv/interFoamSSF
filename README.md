# interFoamSSF
Sharp Surface Force model implemented in interFoam

# Overview
The Sharp Surface Force (SSF) model, proposed by [1], is implemented in interFoam instead of the Continuum Surface Force model as discussed in [2,3]. 

# Description
This version is developed for OpenFOAM 6. Please follow the following step to compile the solver:

Step 1: Go to the directory ```$ cd Libraries/SSF_transportModels/``` and run ``` $ ./Allwclean ``` followed by ``` $ ./Allwmake ```. 

Step 2: Now go the the interFoamSSF folder ```$ cd ../../interFoamSSF/``` and run ``` $ wclean ``` followed by ``` $ wmake ```.

This should have installed interFoamSSF solver in ```FOAM_USER_APPBIN ```.

# Validation 
This implemented Sharp Surface Force model in interFoam, interFoamSSF, is validated for two-dimensional cases of rising bubbles, capillary rise, stationary millimeter and submillimeter sized bubbles. These cases, the solver setting used and the results are described in [2,3].

# Usage
The interFoamSSF solver is released under Creative Commons Attribution-NonCommercial 4.0 International Public License (see LICENSE.txt included included in the main directory). Please feel free to extend, modify and use the solver according to your requirement (only for noncommerical usage), we would only like to ask that you to refer to the papers mentioned in references.

# Related references
1. Raeini, A.Q., Blunt, M.J., Bijeljic, B. 2012. Modelling two-phase flow in porous media at the pore scale using the volume-of-fluid method. Journal of Computational Physics, 231: 5653–5668. [doi:10.1016/j.jcp.2012.04.011](https://doi.org/10.1016/j.jcp.2012.04.011)
2. Vachaparambil, K.J., Einarsrud, K.E. 2019. Comparison of surface tension models for the volume of fluid method. Processes, 7(8): 542. [doi:10.3390/pr7080542](https://doi.org/10.3390/pr7080542)
3. Vachaparambil, K.J., Einarsrud, K.E. 2021. On sharp surface force model: Effect of sharpening coefficient. Experimental and Computational Multiphase Flow, 3: 226–232. [doi:10.1007/s42757-020-0063-5](https://doi.org/10.1007/s42757-020-0063-5)
