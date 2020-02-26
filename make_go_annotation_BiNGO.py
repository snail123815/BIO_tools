# convert go terms got from interproscan to BiNGO compactable files
# https://www.psb.ugent.be/cbd/papers/BiNGO/Customize.html

import numpy as np
import pandas as pd

'''
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
'''

sourceGoFile = '/Users/durand.dc/Documents/works/Resources/Resource_QL109/QL109_GO_20190831.xlsx'
# this file should be dereplicated
source = pd.read_excel(sourceGoFile, header=0, index_col=0)

target = f'{sourceGoFile[:-5]}.tair'
with open(target, 'w') as writer:
    # write header:
    writer.write('(species=Streptomyces roseofaciens)(type=general)(curator=DU)\n')
    # I am not sure but the type seems not important at all, because BiNGO have its own
    # database for the hierarchy. It might be useful only as an identifier if you need to customize
    # the reference to part of the genome/proteome.
    for id, gos in source.iterrows():
        if pd.isna(gos[0]):
            continue
        else:
            gos = gos[0].split('|')
            gos = [go.split(':')[1] for go in gos]
        for go in gos:
            writer.write(f'{id} = {go}\n')
