# To add translation to CDS for the file without it.

from Bio import SeqIO
import sys
import os

seqPath = sys.argv[1]

results = []

for seq in SeqIO.parse(seqPath,'genbank'):
    for feat in seq.features: 
        if feat.type == 'CDS': 
            try: 
                translation = feat.translate(seq,11,cds=True,to_stop=False) 
            except: 
                feat.type = 'CDS_partial' 
                continue 
            feat.qualifiers['translation'] = translation.seq 
    results.append(seq)

outputPath, filename = os.path.split(seqPath)
name, ext = os.path.splitext(filename)
outFilePath = os.path.join(outputPath,f'{name}_filled{ext}')

SeqIO.write(results, outFilePath, 'genbank')