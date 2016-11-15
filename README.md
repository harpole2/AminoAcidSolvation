# AminoAcidSolvation
A set of TCL scripts for calculating water content and adding to bfield of PDB and R scripts for performing analysis

Note: to make this useful for general proteins, the three letter amino acid sequence of specific protein must be added to the 
two auto_read R files

#Pipeline
Start-> waters_near_residue.tcl or coupled_residues_water.tcl > compute_percent_increase_O_C.r > WaterBfield.tcl

