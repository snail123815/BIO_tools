import os

filePath = '/Users/durand.dc/Documents/works/Project_Sporulation/SCO1839_SGR5654/phylogenetic/CDS_DNAseqs_align.txt'
outputPath = f'{os.path.splitext(filePath)[0]}.fasta'
file = open(filePath, 'r')
output = open(outputPath, 'w')

for line in file.readlines():
    if line.startswith('#') or len(line) == 1 or line.startswith('//'):
        continue
    else:
        # print('a')
        # print(line, end='')
        # break
        try:
            name, seq = line.split()
        except:
            print(line)
            break
        # for uniprot alignment:
        # assert name.count('|') == 1
        # name = name.split('|')[1]

        # for DNA alignment there is only one name

        seq = seq.replace('.', '-')
        output.write(f'>{name}\n{seq}\n')

file.close()
output.close()
