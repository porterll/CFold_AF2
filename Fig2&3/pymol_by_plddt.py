import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt



def main():
    data = glob.glob("data.csv")
    df = pd.read_csv(f"{data[0]}")
    df = df.rename(columns=lambda x: x.strip())

    if len(sys.argv) > 1:
        min_plddt = float(sys.argv[1])
    else:
        min_plddt = 0.0
    if len(sys.argv) > 2:
        max_plddt = float(sys.argv[2])
    else:
        max_plddt = 1.0
    

    #don't allow redundant single sequence calculations to appear
    df['single'] = df['file'].str.contains('single', regex=False)
    mask = df['single']
    df = pd.concat([
        df.loc[mask].drop_duplicates(subset = ['single', df.columns[1]], keep='first'),
        df.loc[~mask], ])

    df_bfactors = df[(df["bfactor"] >= min_plddt) & (df["bfactor"] <= max_plddt)]
    df_bfactors = df_bfactors.sort_values(by=['bfactor'], ascending=False)

    files = df_bfactors[df_bfactors.columns[0]].to_list()
    bfactors = df_bfactors["bfactor"].to_list()

    pymol_cmd = "pymol "
    for file in files:
        pymol_cmd += f" {file.strip()}"
    print(pymol_cmd)

if __name__ == "__main__":
    main()
