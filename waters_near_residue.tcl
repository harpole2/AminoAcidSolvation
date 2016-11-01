proc water_around_Manyresidues {args} {
	if {[llength $args] != 4} {
		puts "  Calulates number of water molecules within distance of selection"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res distance"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set dist [lindex $args 3]

        set numFrames [molinfo $simMolID get numframes]

	#loop over all residues selected
	for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
		#atom selection
		set res_tmp [atomselect $simMolID "residue $r and alpha"]


		# get file
		set res_id [$res_tmp get residue]
		set res_name [$res_tmp get resname]
		set file [open "Water-within-$dist-of-residue$res_id-$res_name-notbackbone.dat" w]
		puts $file "frame Water_Molecules Ex_atoms Corrected_water"

		#select water molecules solvating sidechain
		set sel [atomselect $simMolID "same residue as water and within 5 of residue $r and not backbone"] 
	

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
	


	}

}

proc water_around_Manyresidues_allSubunits {args} {
	if {[llength $args] != 6} {
		puts "  Calulates number of water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res distance res_per_subunit _num_subs"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set dist [lindex $args 3]
	set res_per_sub [lindex $args 4]
	set num_subs [lindex $args 5]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {

		#loop over all residues selected
		for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
			#atom selection
			set res_tmp [atomselect $simMolID "residue $r and alpha"]


			# get file
			set res_id [$res_tmp get residue]
			set res_name [$res_tmp get resname]
			set file [open "Water-within-$dist-of-residue$res_id-$res_name-notbackbone.dat" w]
			puts $file "frame Water_Molecules Ex_atoms Corrected_water"

			#select water molecules solvating sidechain
			set sel [atomselect $simMolID "same residue as water and within $dist of residue $r and not backbone"] 
	

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


		}
		set res_first [expr ($res_first + $res_per_sub)]
		set res_last [expr ($res_last + + $res_per_sub)]
	}
}	

proc water_around_Manyresidues_GofRTerminalAtom_allSubunits {args} {
	if {[llength $args] != 5} {
		puts "  Calulates number of water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res res_per_subunit _num_subs"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set res_per_sub [lindex $args 3]
	set num_subs [lindex $args 4]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {

		#loop over all residues selected
		for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
			#atom selection
			set res_tmp [atomselect $simMolID "residue $r and alpha"]
			set res_name [$res_tmp get resname]

			#set distance
			if {$res_name == "ALA"} {
				set dist 5.51
				set atmname CB
			} elseif {$res_name == "ARG"} {
				set dist 5.07
				set atmname CZ
			} elseif {$res_name == "ASN"} {
				set dist 5.71
				set atmname CG
			} elseif {$res_name == "ASP"} {
				set dist 3.15
				set atmname CG
			} elseif {$res_name == "CYS"} {
				set dist 3.9
				set atmname SG
			} elseif {$res_name == "GLN"} {
				set dist 5.8
				set atmname CD
			} elseif {$res_name == "GLU"} {
				set dist 3.15
				set atmname CD
			} elseif {$res_name == "GLY"} {
				set dist 5.37
				set atmname CA
			} elseif {$res_name == "HIS"} {
				set dist 5.55
				set atmname NE2
			} elseif {$res_name == "HSE"} {
				set dist 5.55
				set atmname NE2
			} elseif {$res_name == "HSD"} {
				set dist 5.55
				set atmname NE2
			} elseif {$res_name == "ILE"} {
				set dist 5.55
				set atmname CD
			} elseif {$res_name == "LEU"} {
				set dist 5.48
				set atmname CD2
			} elseif {$res_name == "LYS"} {
				set dist 3.58
				set atmname NZ
			} elseif {$res_name == "MET"} {
				set dist 5.43
				set atmname CE
			} elseif {$res_name == "PHE"} {
				set dist 5.6
				set atmname CZ
			} elseif {$res_name == "PRO"} {
				set dist 5.88
				set atmname CG
			} elseif {$res_name == "SER"} {
				set dist 3.9
				set atmname OG
			} elseif {$res_name == "THR"} {
				set dist 3.85
				set atmname OG1
			} elseif {$res_name == "TRP"} {
				set dist 4.42
				set atmname NE1
			} elseif {$res_name == "TYR"} {
				set dist 4.08
				set atmname OH
			} elseif {$res_name == "VAL"} {
				set dist 5.28
				set atmname CG1
			}

			# get file
			set res_id [$res_tmp get residue]
			set res_name [$res_tmp get resname]
			set file [open "Water-within-GofR-of-residue$res_id-$res_name-notbackbone.dat" w]
			puts $file "frame Water_Molecules Ex_atoms Corrected_water"

			#select water molecules solvating sidechain
			set sel [atomselect $simMolID "same residue as water and within $dist of residue $r and name $atmname"] 
	

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


		}
		set res_first [expr ($res_first + $res_per_sub)]
		set res_last [expr ($res_last + + $res_per_sub)]
	}
}

proc water_around_Manyresidues_GofR_allSubunits {args} {
	if {[llength $args] != 5} {
		puts "  Calulates number of water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res res_per_subunit _num_subs"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set res_per_sub [lindex $args 3]
	set num_subs [lindex $args 4]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {

		#loop over all residues selected
		for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
			#atom selection
			set res_tmp [atomselect $simMolID "residue $r and alpha"]
			set res_name [$res_tmp get resname]

			#set distance
			if {$res_name == "ALA"} {
				set dist 5.51

			} elseif {$res_name == "ARG"} {
				set dist 5.07

			} elseif {$res_name == "ASN"} {
				set dist 5.71

			} elseif {$res_name == "ASP"} {
				set dist 3.15

			} elseif {$res_name == "CYS"} {
				set dist 3.9

			} elseif {$res_name == "GLN"} {
				set dist 5.8

			} elseif {$res_name == "GLU"} {
				set dist 3.15

			} elseif {$res_name == "GLY"} {
				set dist 5.37

			} elseif {$res_name == "HIS"} {
				set dist 5.55

			} elseif {$res_name == "HSD"} {
				set dist 5.55

			} elseif {$res_name == "HSE"} {
				set dist 5.55

			} elseif {$res_name == "ILE"} {
				set dist 5.55

			} elseif {$res_name == "LEU"} {
				set dist 5.48

			} elseif {$res_name == "LYS"} {
				set dist 3.58

			} elseif {$res_name == "MET"} {
				set dist 5.43

			} elseif {$res_name == "PHE"} {
				set dist 5.6

			} elseif {$res_name == "PRO"} {
				set dist 5.88

			} elseif {$res_name == "SER"} {
				set dist 3.9

			} elseif {$res_name == "THR"} {
				set dist 3.85

			} elseif {$res_name == "TRP"} {
				set dist 4.42

			} elseif {$res_name == "TYR"} {
				set dist 4.08

			} elseif {$res_name == "VAL"} {
				set dist 5.28

			}

			# get file
			set res_id [$res_tmp get residue]
			set res_name [$res_tmp get resname]
			set file [open "Water-within-GofR-of-residue$res_id-$res_name-notbackbone.dat" w]
			puts $file "frame Water_Molecules Ex_atoms Corrected_water"

			#select water molecules solvating sidechain
			set sel [atomselect $simMolID "same residue as water and within $dist of residue $r and not backbone"] 
	

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


		}
		set res_first [expr ($res_first + $res_per_sub)]
		set res_last [expr ($res_last + + $res_per_sub)]
	}
}

proc water_around_Manyresidues_GofR_allSubunits_newR {args} {
	if {[llength $args] != 5} {
		puts "  Calulates number of water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res res_per_subunit _num_subs"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set res_per_sub [lindex $args 3]
	set num_subs [lindex $args 4]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {

		#loop over all residues selected
		for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
			#atom selection
			set res_tmp [atomselect $simMolID "residue $r and alpha"]
			set res_name [$res_tmp get resname]

			#set distance
			if {$res_name == "ALA"} {
				set dist 5.45

			} elseif {$res_name == "ARG"} {
				set dist 4.75

			} elseif {$res_name == "ASN"} {
				set dist 5.49

			} elseif {$res_name == "ASP"} {
				set dist 3.19

			} elseif {$res_name == "CYS"} {
				set dist 3.5

			} elseif {$res_name == "GLN"} {
				set dist 5.38

			} elseif {$res_name == "GLU"} {
				set dist 3.18

			} elseif {$res_name == "GLY"} {
				set dist 5.49

			} elseif {$res_name == "HIS"} {
				set dist 5.48

			} elseif {$res_name == "HSE"} {
				set dist 5.48

			} elseif {$res_name == "HSD"} {
				set dist 5.48

			} elseif {$res_name == "ILE"} {
				set dist 5.45

			} elseif {$res_name == "LEU"} {
				set dist 5.42

			} elseif {$res_name == "LYS"} {
				set dist 3.58

			} elseif {$res_name == "MET"} {
				set dist 5.18

			} elseif {$res_name == "PHE"} {
				set dist 5.28

			} elseif {$res_name == "PRO"} {
				set dist 5.65

			} elseif {$res_name == "SER"} {
				set dist 3.5

			} elseif {$res_name == "THR"} {
				set dist 3.48

			} elseif {$res_name == "TRP"} {
				set dist 3.75

			} elseif {$res_name == "TYR"} {
				set dist 3.68

			} elseif {$res_name == "VAL"} {
				set dist 5.35

			}

			# get file
			set res_id [$res_tmp get residue]
			set res_name [$res_tmp get resname]
			set file [open "Water-within-GofR-of-residue$res_id-$res_name-notbackbone.dat" w]
			puts $file "frame Water_Molecules Ex_atoms Corrected_water"

			#select water molecules solvating sidechain
			set sel [atomselect $simMolID "same residue as water and within $dist of residue $r and not backbone"] 
	

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


		}
		set res_first [expr ($res_first + $res_per_sub)]
		set res_last [expr ($res_last + + $res_per_sub)]
	}
}

proc water_around_Manyresidues_GofRTerminalAtom_allSubunits_newR {args} {
	if {[llength $args] != 5} {
		puts "  Calulates number of water molecules within distance of selection all subunits"
        	puts "  Usage: water_around_residue sim_mol_id first_res last_res res_per_subunit _num_subs"
		puts "  Note: do not use residue first_res, just the number"
		puts "  Note: this will go through all residues between first_res last_res"
        	error ""
	}
	# Parse the arguments.
	set simMolID [lindex $args 0]
	set res_first [lindex $args 1]
	set res_last [lindex $args 2]
	set res_per_sub [lindex $args 3]
	set num_subs [lindex $args 4]

        set numFrames [molinfo $simMolID get numframes]


	#loop over all four subunits
	for {set c 1} {$c <= $num_subs} {incr c} {

		#loop over all residues selected
		for {set r $res_first} {$r <= $res_last} {incr r } {
	
		
			#atom selection
			set res_tmp [atomselect $simMolID "residue $r and alpha"]
			set res_name [$res_tmp get resname]

			#set distance
			if {$res_name == "ALA"} {
				set dist 5.45
				set atmname CB
			} elseif {$res_name == "ARG"} {
				set dist 4.75
				set atmname CZ
			} elseif {$res_name == "ASN"} {
				set dist 5.49
				set atmname CG
			} elseif {$res_name == "ASP"} {
				set dist 3.19
				set atmname CG
			} elseif {$res_name == "CYS"} {
				set dist 3.5
				set atmname SG
			} elseif {$res_name == "GLN"} {
				set dist 5.38
				set atmname CD
			} elseif {$res_name == "GLU"} {
				set dist 3.18
				set atmname CD
			} elseif {$res_name == "GLY"} {
				set dist 5.49
				set atmname CA
			} elseif {$res_name == "HIS"} {
				set dist 5.48
				set atmname NE2
			} elseif {$res_name == "HSE"} {
				set dist 5.48
				set atmname NE2
			} elseif {$res_name == "HSD"} {
				set dist 5.48
				set atmname NE2
			} elseif {$res_name == "ILE"} {
				set dist 5.45
				set atmname CD
			} elseif {$res_name == "LEU"} {
				set dist 5.42
				set atmname CD2
			} elseif {$res_name == "LYS"} {
				set dist 3.58
				set atmname NZ
			} elseif {$res_name == "MET"} {
				set dist 5.18
				set atmname CE
			} elseif {$res_name == "PHE"} {
				set dist 5.28
				set atmname CZ
			} elseif {$res_name == "PRO"} {
				set dist 5.65
				set atmname CG
			} elseif {$res_name == "SER"} {
				set dist 3.5
				set atmname OG
			} elseif {$res_name == "THR"} {
				set dist 3.48
				set atmname OG1
			} elseif {$res_name == "TRP"} {
				set dist 3.75
				set atmname NE1
			} elseif {$res_name == "TYR"} {
				set dist 3.68
				set atmname OH
			} elseif {$res_name == "VAL"} {
				set dist 5.35
				set atmname CG1
			}

			# get file
			set res_id [$res_tmp get residue]
			set res_name [$res_tmp get resname]
			set file [open "Water-within-GofR-of-residue$res_id-$res_name-notbackbone.dat" w]
			puts $file "frame Water_Molecules Ex_atoms Corrected_water"

			#select water molecules solvating sidechain
			set sel [atomselect $simMolID "same residue as water and within $dist of residue $r and name $atmname"] 
	

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


		}
		set res_first [expr ($res_first + $res_per_sub)]
		set res_last [expr ($res_last + + $res_per_sub)]
	}
}

