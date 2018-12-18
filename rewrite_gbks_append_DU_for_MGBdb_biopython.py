import os
import shutil
from Bio import SeqIO
# inspired by Alexander

def rewrite_gbk_append(infolder,outfolder):
    genomes = [file for file in os.listdir(infolder) if os.path.splitext(file)[1] in ['.gb','.gbk','.genbank'] ]
    geneSeen = {}
    proteinSeen = {}
    totalAll = 0
    nonUniqueAll = 0
    for genome in genomes:
        sourceSeqs = SeqIO.parse(os.path.join(infolder,genome), 'genbank')
        newSeqs = []
        total = 0
        nonUnique = 0
        for seq in sourceSeqs:
            for feat in seq.features:
                if feat.type == 'CDS':
                    total += 1
                    if 'gene' in feat.qualifiers:
                        oriGeneName = feat.qualifiers['gene'][0].replace('.', '_')
                        feat.qualifiers['gene'][0] = oriGeneName
                        if oriGeneName not in geneSeen:
                            geneSeen[oriGeneName] = 0
                            feat.qualifiers['gene'][0] = oriGeneName
                        else:
                            geneSeen[oriGeneName] += 1
                            newGeneName = '_'.join([oriGeneName, str(geneSeen[oriGeneName])])
                            feat.qualifiers['gene'][0] = newGeneName
                            nonUnique += 1
                    if 'protein_id' in feat.qualifiers:
                        oriProteinName = feat.qualifiers['protein_id'][0].replace('.', '_')
                        if oriProteinName not in proteinSeen:
                            proteinSeen[oriProteinName] = 0
                            feat.qualifiers['protein_id'] = oriProteinName
                        else:
                            proteinSeen[oriProteinName] += 1
                            newProteinName = '_'.join([oriProteinName,str(proteinSeen[oriProteinName])])
                            feat.qualifiers['protein_id'] = newProteinName
                            nonUnique += 1
            newSeqs.append(seq)
        print('In %s found %i accessions, %i of which were Non-unique' %(genome,total,nonUnique))
        print('Writing converted file...')
        SeqIO.write(newSeqs, os.path.join(outfolder,genome), 'genbank')
        totalAll += total
        nonUniqueAll += nonUnique
    return totalAll, nonUniqueAll


if __name__ == '__main__':
    import sys
    args = sys.argv
    if len(args) < 3:
        print('\nUSAGE: python rewrite_gbks.py infolder outfolder\n')
        exit()
    infolder = args[1]
    outfolder = args[2]
    if os.path.isdir(outfolder):
        print('\nWarning: overwriting outputfolder. Continue?\n')
        i = input('y/n:')
        if i.lower() != 'y':
            exit()
        shutil.rmtree(outfolder)
    os.mkdir(outfolder)
    totalAll, nonUniqueAll = rewrite_gbk_append(infolder,outfolder)
    print('\nFinished converting.\nFound %i accessions, %i of which were Non-unique.\n' %(totalAll,nonUniqueAll))
