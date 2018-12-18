import csv
import os
dir_path = os.path.dirname(os.path.realpath(__file__))

inputFile = os.path.join(dir_path,'/Users/durand.dc/Desktop/ChIP1839/temp/fimo.csv')
outputFile = os.path.join(dir_path, '/Users/durand.dc/Desktop/ChIP1839/temp/fimo.fasta')

with open(inputFile, 'r') as inputHandle:
    reader = csv.reader(inputHandle)
    with open(outputFile, 'w') as outputHandle:
        for row in reader:
            outputHandle.write('>{0}\n{1}\n'.format(row[0], row[1]))
