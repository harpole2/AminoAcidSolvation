proc water_around_Coupledresidues_GofRTerminalAtom_allSubunits {args} {
	if {[llength $args] != 5} {
		puts "  Calulates number of unique water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res second_res res_per_subunit _num_subs"
		puts "  Note: this is intended to work for only 2 neighboring residues"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_one [lindex $args 1]
	set res_two [lindex $args 2]
	set res_per_sub [lindex $args 3]
	set num_subs [lindex $args 4]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {
	
		
		#atom selection
		set res_tmp_1 [atomselect $simMolID "residue $res_one and alpha"]
		set resname_1 [$res_tmp_1 get resname]

		set res_tmp_2 [atomselect $simMolID "residue $res_two and alpha"]
		set resname_2 [$res_tmp_2 get resname]

		#set distance1
		if {$resname_1 == "ALA"} {
			set dist1 5.51
			set atmname1 CB
		} elseif {$resname_1 == "ARG"} {
			set dist1 5.07
			set atmname1 CZ
		} elseif {$resname_1 == "ASN"} {
			set dist1 5.71
			set atmname1 CG
		} elseif {$resname_1 == "ASP"} {
			set dist1 3.15
			set atmname1 CG
		} elseif {$resname_1 == "CYS"} {
			set dist1 3.9
			set atmname1 SG
		} elseif {$resname_1 == "GLN"} {
			set dist1 5.8
			set atmname1 CD
		} elseif {$resname_1 == "GLU"} {
			set dist1 3.15
			set atmname1 CD
		} elseif {$resname_1 == "GLY"} {
			set dist1 5.37
			set atmname1 CA
		} elseif {$resname_1 == "HIS"} {
			set dist1 5.55
			set atmname1 NE2
		} elseif {$resname_1 == "HSE"} {
			set dist1 5.55
			set atmname1 NE2
		} elseif {$resname_1 == "HSD"} {
			set dist1 5.55
			set atmname1 NE2
		} elseif {$resname_1 == "ILE"} {
			set dist1 5.55
			set atmname1 CD
		} elseif {$resname_1 == "LEU"} {
			set dist1 5.48
			set atmname1 CD2
		} elseif {$resname_1 == "LYS"} {
			set dist1 3.58
			set atmname1 NZ
		} elseif {$resname_1 == "MET"} {
			set dist1 5.43
			set atmname1 CE
		} elseif {$resname_1 == "PHE"} {
			set dist1 5.6
			set atmname1 CZ
		} elseif {$resname_1 == "PRO"} {
			set dist1 5.88
			set atmname1 CG
		} elseif {$resname_1 == "SER"} {
			set dist1 3.9
			set atmname1 OG
		} elseif {$resname_1 == "THR"} {
			set dist1 3.85
			set atmname1 OG1
		} elseif {$resname_1 == "TRP"} {
			set dist1 4.42
			set atmname1 NE1
		} elseif {$resname_1 == "TYR"} {
			set dist1 4.08
			set atmname1 OH
		} elseif {$resname_1 == "VAL"} {
			set dist1 5.28
			set atmname1 CG1
		}

		#set distance2
		if {$resname_2 == "ALA"} {
			set dist2 5.51
			set atmname2 CB
		} elseif {$resname_2 == "ARG"} {
			set dist2 5.07
			set atmname2 CZ
		} elseif {$resname_2 == "ASN"} {
			set dist2 5.71
			set atmname2 CG
		} elseif {$resname_2 == "ASP"} {
			set dist2 3.15
			set atmname2 CG
		} elseif {$resname_2 == "CYS"} {
			set dist2 3.9
			set atmname2 SG
		} elseif {$resname_2 == "GLN"} {
			set dist2 5.8
			set atmname2 CD
		} elseif {$resname_2 == "GLU"} {
			set dist2 3.15
			set atmname2 CD
		} elseif {$resname_2 == "GLY"} {
			set dist2 5.37
			set atmname2 CA
		} elseif {$resname_2 == "HIS"} {
			set dist2 5.55
			set atmname2 NE2
		} elseif {$resname_2 == "HSE"} {
			set dist2 5.55
			set atmname2 NE2
		} elseif {$resname_2 == "HSD"} {
			set dist2 5.55
			set atmname2 NE2
		} elseif {$resname_2 == "ILE"} {
			set dist2 5.55
			set atmname2 CD
		} elseif {$resname_2 == "LEU"} {
			set dist2 5.48
			set atmname2 CD2
		} elseif {$resname_2 == "LYS"} {
			set dist2 3.58
			set atmname2 NZ
		} elseif {$resname_2 == "MET"} {
			set dist2 5.43
			set atmname2 CE
		} elseif {$resname_2 == "PHE"} {
			set dist2 5.6
			set atmname2 CZ
		} elseif {$resname_2 == "PRO"} {
			set dist2 5.88
			set atmname2 CG
		} elseif {$resname_2 == "SER"} {
			set dist2 3.9
			set atmname2 OG
		} elseif {$resname_2 == "THR"} {
			set dist2 3.85
			set atmname2 OG1
		} elseif {$resname_2 == "TRP"} {
			set dist2 4.42
			set atmname2 NE1
		} elseif {$resname_2 == "TYR"} {
			set dist2 4.08
			set atmname2 OH
		} elseif {$resname_2 == "VAL"} {
			set dist2 5.28
			set atmname2 CG1
		}

		# get file
		set res_id1 [$res_tmp_1 get residue]
		set res_name1 [$res_tmp_1 get resname]
		set res_id2 [$res_tmp_2 get residue]
		set res_name2 [$res_tmp_2 get resname]
		set file [open "Water-within-GofR-of-residue$res_id1-$res_name1-or-residue$res_id2-$res_name2-Coupled.dat" w]
		puts $file "frame Water_Molecules Ex_atoms Corrected_water"

		#select water molecules solvating sidechain
		set sel [atomselect $simMolID "(same residue as water and within $dist1 of residue $res_one and name $atmname1) or (same residue as water and within $dist2 of residue $res_two and name $atmname2)"] 
	

		# Selection loop
		for {set i 0 } {$i < $numFrames } { incr i } {
			$sel frame $i
			$sel update
			set number [$sel num]
			#num gives atoms, we want residues
			set numberW [expr ($number/3)]
			set numberR [expr ($number%3)]
			# start corrected number with original number
			set numberWp $numberW
			#if there are extra atoms because a hydrogen wasn't in same residue and therefore water wasn't counted
			#add one more water molecule
			if {$numberR > 0} {
				set numberWp [expr $numberW +1]
			}
			puts $file "$i $numberW $numberR $numberWp"
		}
		close $file
		$sel delete
		unset sel


		
		set res_one [expr ($res_one + $res_per_sub)]
		set res_two [expr ($res_two + + $res_per_sub)]
	}
}


proc water_around_Coupledresidues_dist_allSubunits {args} {
	if {[llength $args] != 6} {
		puts "  Calulates number of unique water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res second_res  distance res_per_subunit _num_subs"
		puts "  Note: this is intended to work for only 2 neighboring residues"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_one [lindex $args 1]
	set res_two [lindex $args 2]
	set dist [lindex $args 3]
	set res_per_sub [lindex $args 4]
	set num_subs [lindex $args 5]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {
	
		
		#atom selection
		set res_tmp_1 [atomselect $simMolID "residue $res_one and alpha"]
		set resname_1 [$res_tmp_1 get resname]

		set res_tmp_2 [atomselect $simMolID "residue $res_two and alpha"]
		set resname_2 [$res_tmp_2 get resname]

		#set distance1

		# get file
		set res_id1 [$res_tmp_1 get residue]
		set res_name1 [$res_tmp_1 get resname]
		set res_id2 [$res_tmp_2 get residue]
		set res_name2 [$res_tmp_2 get resname]
		set file [open "Water-within-GofR-of-residue$res_id1-$res_name1-or-residue$res_id2-$res_name2-Coupled.dat" w]
		puts $file "frame Water_Molecules Ex_atoms Corrected_water"

		#select water molecules solvating sidechain
		set sel [atomselect $simMolID "(same residue as water and within $dist of residue $res_one and not backbone) or (same residue as water and within $dist of residue $res_two and not backbone)"] 
	

		# Selection loop
		for {set i 0 } {$i < $numFrames } { incr i } {
			$sel frame $i
			$sel update
			set number [$sel num]
			#num gives atoms, we want residues
			set numberW [expr ($number/3)]
			set numberR [expr ($number%3)]
			# start corrected number with original number
			set numberWp $numberW
			#if there are extra atoms because a hydrogen wasn't in same residue and therefore water wasn't counted
			#add one more water molecule
			if {$numberR > 0} {
				set numberWp [expr $numberW +1]
			}
			puts $file "$i $numberW $numberR $numberWp"
		}
		close $file
		$sel delete
		unset sel


		
		set res_one [expr ($res_one + $res_per_sub)]
		set res_two [expr ($res_two + + $res_per_sub)]
	}
}

