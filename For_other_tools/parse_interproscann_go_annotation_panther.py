import pandas as pd
import argparse
import os


desctiptionStr = '''
Extract and de-replicate PANTHER ID for PANTHER analysis (Overrepresentation)

Output:

col 1: sequence ID
col 2: PANTHER accession (PTHRnnnnn for family HMMs, PTHRnnnnn:SFnn for subfamilies)


'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract and de-replicate GO terms for BiNGO reference')
    parser.add_argument('path', help="path to interproscan output. needs '-appl PANTHER'")
    parser.add_argument('--species', help='species', default='Not_set')
    parser.add_argument(
        '--curator', help='curator, generator of the annotation', default='Not_set')
    args = parser.parse_args()
    interproscanOutTsv = args.path.strip()
    species = args.species
    curator = args.curator
    convertedFile = f'{os.path.splitext(interproscanOutTsv)[0]}_pantherRef.tsv'
    source = pd.read_csv(interproscanOutTsv, sep='\t', header=None, usecols=range(9), index_col=None)
    # print(source.head(10))
    corrdict = {} # id:panacc
    for i,r in source.iterrows():
        proid = r[0]
        if proid not in corrdict:
            corrdict[proid] = ["-", "-"] # retrieve none mapped for further notice
        if r[3] != "PANTHER":
            continue
        panacc = r[4]
        alrange = f"{r[6]}-{r[7]}"
        if proid in corrdict and corrdict[proid][0] != "-":
            if ':' in panacc:
                if ':' in corrdict[proid][0] and panacc != corrdict[proid][0]:
                    print('EXCEPTION:')
                    print(f'{proid}: {corrdict[proid][0]}[{corrdict[proid][1]}], {panacc}[{alrange}]')
                    exit()
                else:
                    corrdict[proid] = [panacc, alrange]
            elif ':' not in corrdict[proid][0] and panacc != corrdict[proid][0]:
                print('EXCEPTION:')
                print(f'{proid}: {corrdict[proid][0]}[{corrdict[proid][1]}], {panacc}[{alrange}]')
                exit()
        else:
            corrdict[proid] = [panacc, alrange]
    with open(convertedFile, 'w') as writer:
        for proid in corrdict:
            if corrdict[proid][0] != '-': # PANTHER does not accept none mapped
                writer.write("\t".join([proid,corrdict[proid][0]]))
                writer.write('\n')

    print(convertedFile)
