###############################################
# Scripted by Chao DU                         #
# c.du@biology.leidenuniv.nl                  #
# durand.dc@hotmail.com                       #
###############################################

import os
import shutil
from Bio import SeqIO
import threading


def rewrite_gbk_append(genome, infolder, counter):
    sourceSeqs = SeqIO.parse(os.path.join(infolder, genome), 'genbank')
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
                        newGeneName = '_'.join(
                            [oriGeneName, str(geneSeen[oriGeneName])])
                        feat.qualifiers['gene'][0] = newGeneName
                        nonUnique += 1
                if 'protein_id' in feat.qualifiers:
                    oriProteinName = feat.qualifiers['protein_id'][0].replace(
                        '.', '_')
                    if oriProteinName not in proteinSeen:
                        proteinSeen[oriProteinName] = 0
                        feat.qualifiers['protein_id'] = oriProteinName
                    else:
                        proteinSeen[oriProteinName] += 1
                        newProteinName = '_'.join(
                            [oriProteinName, str(proteinSeen[oriProteinName])])
                        feat.qualifiers['protein_id'] = newProteinName
                        nonUnique += 1
        newSeqs.append(seq)
    print('In %s found %i accessions, %i of which were Non-unique' %
          (genome, total, nonUnique))
    counter[genome] = [total, nonUnique, newSeqs]
# rewrite_gbk_append


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
    geneSeen = {}
    proteinSeen = {}
    counter = {}
    totalAll = 0
    nonUniqueAll = 0
    genomes = [file for file in os.listdir(infolder) if
               os.path.splitext(file)[1] in ['.gb', '.gbk', '.genbank']]
    threads = {}
    for genome in genomes:
        threads[genome] = threading.Thread(
            target=rewrite_gbk_append, args=(genome, infolder, counter))
        threads[genome].start()
    for genome in genomes:
        threads[genome].join()
    for genome in genomes:
        total, nonUnique, newSeqs = counter[genome]
        outfile = os.path.join(outfolder, genome)
        totalAll += total
        nonUniqueAll += nonUnique
        print('Writing %s...' % (genome))
        SeqIO.write(newSeqs, outfile, 'genbank')
    print('\nFinished converting.\nFound %i accessions, %i of which were Non-unique.\n' %
          (totalAll, nonUniqueAll))
