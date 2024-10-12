#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from tmalign import save_pdb_without_ter,alignwithanymethod,tmalign,tmscore,mmalign
import pymol
import sys
import glob
import pandas as pd
import numpy as np
from collections import defaultdict

def CL_input():
    """
    Parse command line arguments that are being passed in
    """

    error_message = """
Missing command line arguments!')
Available flags:
    -path  ####  |  Directory containing CF-random runs
    -ref_1 ####  |  first  reference model to compare all structures
    -ref_2 ####  |  second reference model to compare all structures
    """
    if len(sys.argv) < 5:
        print(error_message)
        sys.exit()

    options = ['-path', '-ref_1', '-ref_2']
    for i in sys.argv[1:]:
        if '-' in i and i not in options:
            print(f"This is not an option: {i}")
            print(error_message)
            sys.exit()
    path =  sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-path' == _][0] + 1]
    ref_1 = sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-ref_1' == _][0] + 1]
    ref_2 = sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-ref_2' == _][0] + 1]

    return path, ref_1, ref_2


def main():

    path, ref_1, ref_2 = CL_input()
    #path should be the outermost directory containing pdb files
    #-path .  <-- get every structure in every subdirectory
    pdb_structures = sorted(glob.glob('**/*.pdb', recursive=True)) #<- careful here this will grab everything
    pdb_structures = [structure for structure in pdb_structures if not (re.search(fr'({structure})', ref_2) or re.search(fr'({structure})', ref_1))]

    #upload model and template pdbs to pymol session
    reference_1 = ref_1
    pymol.cmd.load(reference_1)

    reference_2 = ref_2
    pymol.cmd.load(reference_2)

    data = defaultdict(list)    
    for pdb in pdb_structures:
        #retrieve names from pymol session
        pymol.cmd.load(pdb)
        data["file"].append(pdb)
        names = pymol.cmd.get_names()
        print(names)
        # calculate tmscore
        t1_sub = tmscore(f"{names[0]} & i. 118-155",f"{names[2]} & i. 118-155", quiet=1, exe='/Users/schaferjw/pymol_setup/TMalign')
        data[f"{names[0]}_tm_sub"].append(t1_sub)
        t1 = tmscore(names[0], names[2], quiet=1, exe='/Users/schaferjw/pymol_setup/TMalign')
        data[f"{names[0]}_tm"].append(t1)

        t2_sub = tmscore(f"{names[1]} & i. 118-155",f"{names[2]} & i. 118-155", quiet=1, exe='/Users/schaferjw/pymol_setup/TMalign')
        data[f"{names[1]}_tm_sub"].append(t2_sub)
        t2 = tmscore(names[1], names[2], quiet=1, exe='/Users/schaferjw/pymol_setup/TMalign')
        data[f"{names[1]}_tm"].append(t2)


        #calculate rmsd
        rmsd1 = pymol.cmd.align(f"{names[0]} & i. 118-155",f"{names[2]} & i. 118-155")
        data[f"{names[0]}_rmsd"].append(rmsd1[0])
        data[f"{names[0]}_rmsd_res_aligned"].append(rmsd1[-1])
        data[f"{names[0]}_rmsd_before_align"].append(rmsd1[3])


        rmsd2 = pymol.cmd.align(f"{names[1]} & i. 118-155",f"{names[2]} & i. 118-155")
        data[f"{names[1]}_rmsd"].append(rmsd2[0])
        data[f"{names[1]}_rmsd_res_aligned"].append(rmsd2[-1])
        data[f"{names[1]}_rmsd_before_align"].append(rmsd2[3])

        #grab plddt values
        predicted_model = pymol.cmd.get_model(f"{names[2]} and name CA")
        bf = [atom.b for atom in predicted_model.atom if atom.name == "CA"] #<redundant but safe
        data["bfactor"].append(np.average(np.array(bf)))

        #close current pymol session and remove all data from memory
        pymol.cmd.delete(names[2])

    df = pd.DataFrame.from_dict(data)
    df.to_csv("data.csv", index=False)


if __name__ == "__main__":
    """
        pymol.cmd.align returns a list with 7 items:

        RMSD after refinement
        Number of aligned atoms after refinement
        Number of refinement cycles
        RMSD before refinement
        Number of aligned atoms before refinement
        Raw alignment score
        Number of residues aligned
    """

    main()
