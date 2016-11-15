#Use b_field table created  in compute_percent_increase_O_C.r

proc BfieldWaterColor {args} {
	if {[llength $args] != 2} {
		puts "  Puts water content in Bfield of PDB for visualization"
        	puts "  Usage: BfieldWaterColor sim_mol_id residue_bfield_file"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set bfile [lindex $args 1]

	#read file
	set inFile [open $bfile r]
    	set numAtoms 0
    	while { [gets $inFile line] >= 0 } {

	        set splitLine [split [string trim $line]]
       	 	set resNum($numAtoms) [lindex $splitLine 0]
        	set bfield($numAtoms) [lindex $splitLine end]

		incr numAtoms
	}
	close $inFile

    	set sel [atomselect $simMolID "all"]
    	$sel set beta 0.0
    	$sel delete
    	for {set i 0} {$i < $numAtoms} {incr i} {
        	set sel [atomselect $simMolID "residue $resNum($i)"]
        	$sel set beta $bfield($i)
        	$sel delete
	}
}

