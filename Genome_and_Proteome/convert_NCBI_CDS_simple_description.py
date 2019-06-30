# This script is to generate a simplified fasta header. The simplified fasta file can be used in PLGS database as 'long description' header database


from collections import Counter

file = "/Users/durand.dc/Downloads/sequence.txt"
target = "/Users/durand.dc/Downloads/SCP1_prot.fasta"

IDList = []

with open(file, 'r') as input:
    with open(target, 'w') as output:
        for l in input.readlines():
            if l.startswith('>'):
                tagList = l.split('[')
                for tag in tagList:
                    if tag.startswith('locus_tag='):
                        protID = tag.split('=')[1][:-2]
                        IDList.append(protID)
                        n = Counter(IDList)[protID]
                        if n > 1:
                            protID = f'{protID}_{n-1}'
                    elif tag.startswith('protein='):
                        protName = tag.split('=')[1][:-2]
                newTag = f'>{protID} {protName}\n'
                output.write(newTag)
            else:
                output.write(l)

duplicateIDs = [item for item, count in Counter(IDList).items() if count > 1]
if len(duplicateIDs) == 0:
    print('No duplicate IDs.')
else:
    print('Duplicates are:')
    print(duplicateIDs)
    print('Numbers appended to the end of these duplicate IDs.')
