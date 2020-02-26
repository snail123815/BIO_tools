import pandas as pd
import argparse
import os


desctiptionStr = '''
Extract and de-replicate GO terms for BiNGO reference:
The format of custom annotation/ontology files in BiNGO is the same as the Cytoscape annotation and ontology file formats. To make a custom annotation file, just parse your annotation into the following form :

(species=Saccharomyces cerevisiae)(type=Biological Process)(curator=GO)
YAL001C = 0006384
YAL002W = 0045324
YAL002W = 0045324
YAL003W = 0006414
YAL004W = 0000004
YAL005C = 0006616
YAL005C = 0006457
YAL005C = 0000060
YAL007C = 0006888
YAL008W = 0000004
...
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract and de-replicate GO terms for BiNGO reference')
    parser.add_argument('path')
    parser.add_argument('--species', help='species', default='Not_set')
    parser.add_argument(
        '--curator', help='curator, generator of the annotation', default='Not_set')
    args = parser.parse_args()
    interproscanOutTsv = args.path.strip()
    species = args.species
    curator = args.curator
    convertedFile = f'{os.path.splitext(interproscanOutTsv)[0]}_GO.tair'
    source = pd.read_csv(interproscanOutTsv, sep='\t', header=None, usecols=[0, 13])
    # print(source.head(10))
    with open(convertedFile, 'w') as writer:
        # write header:
        writer.write(f'(species={species})(type=general)(curator={curator})\n')
        # I am not sure but the type seems not important at all, because BiNGO have its own
        # database for the hierarchy. It might be useful only as an identifier if you need to customize
        # the reference to part of the genome/proteome.
        for num, line in source.iterrows():
            if pd.isna(line.loc[13]):
                continue
            else:
                gos = line.loc[13].split('|')
                gos = [go.split(':')[1] for go in gos]
                id = line.loc[0]
            for go in gos:
                writer.write(f'{id} = {go}\n')

    print(convertedFile)
