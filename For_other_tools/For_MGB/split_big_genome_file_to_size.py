import os
import shutil
from Bio import SeqIO


def shrinkFileSize(infolder, outfolder, sizeLimit=100000000):
    genomes = [file for file in os.listdir(infolder) if
               os.path.splitext(file)[1] in ['.gb', '.gbk', '.genbank']]
    for genome in genomes:
        genomeFile = os.path.join(infolder, genome)
        fileSize = os.path.getsize(genomeFile)
        if fileSize > sizeLimit:
            print(f'Processing {genome}')
            backupFile = os.path.join(outfolder, genome)
            shutil.move(genomeFile, backupFile)
            bigSeqs = SeqIO.parse(backupFile, 'genbank')
            length = 0
            count = 0
            newSeqs = []

            for seq in bigSeqs:
                length += len(seq)
                newSeqs.append(seq)
                if length < sizeLimit * 0.4 / 2:
                    pass
                else:
                    outputFile = os.path.join(infolder, "%s_%i%s" % (
                        os.path.splitext(genome)[0], count, os.path.splitext(genome)[1]))
                    print(f'Reach limit, length = {length}')
                    SeqIO.write(newSeqs, outputFile, 'genbank')
                    newSeqs = []
                    count += 1
                    length = 0
            print(f'Reach end, length = {length}\n')
            if length != 0:
                outputFile = os.path.join(infolder, "%s_%i%s" % (
                    os.path.splitext(genome)[0], count, os.path.splitext(genome)[1]))
                SeqIO.write(newSeqs, outputFile, 'genbank')


# shrinkFileSize
if __name__ == '__main__':
    import sys
    args = sys.argv
    if len(args) < 3:
        print('\nUSAGE: python shrinked.py infolder outfolder\n')
        exit()
    infolder = args[1]
    outfolder = args[2]
    if os.path.isdir(outfolder):
        print('\nWarning: overwriting backup. Continue?\n')
        i = input('y/n:')
        if i.lower() != 'y':
            exit()
        shutil.rmtree(outfolder)
    os.mkdir(outfolder)
    shrinkFileSize(infolder, outfolder)
