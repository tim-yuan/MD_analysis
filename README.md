# MD_analysis

All the codes are written by Tim Yuan

For academic purpose only

Use it at your own risk

## density of z
This code calculates the density of a selected group in the z direction

Required: GROMACS 2018 or above, a trajectory file (xtc or trr), a gro/tpr file

## q3
This code is written based on Li T., et al, PCCP, 13 44 19807--19813, 2011.

This code is based on GROMACS 2018 or more recent versions, a make file is required from GROMACS template.

It takes the index of OW atoms, and identify the solid-like clusters using Steinhardt order parameter.

Note that, the code has a section to compute the histogram of the q3 distribution (uncomment the necessary lines). Please set the cutoff distance to 0.35nm in order to have an exact match to the histogram published by Li, et al, 2011. 


Required: GROMACS 2018 or above, a trajectory file (xtc or trr), a gro/tpr file

## instantaneous liquid interface
This code is written by Tim Yuan, based on Willard A. P. and Chandler D. JPCB, 114,5,1954-1958,2010


