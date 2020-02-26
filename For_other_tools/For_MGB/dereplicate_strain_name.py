###############################################
# Scripted by Chao DU                         #
# c.du@biology.leidenuniv.nl                  #
# durand.dc@hotmail.com                       #
###############################################


# Open log file created when making multi gene blast MGB database,
# count how many strains there are, and sort base on type -
# WGS, plasmid, shotgun plasmid, chromosome, cds

import pandas as pd


table = pd.read_csv('/Users/durand.dc/Documents/works/Misc.files/Streptomyces_genomes/MGB_DB/all_actino/actino_gbk_all_descrs.txt',
                    sep='\t', header=None, index_col=0)

names = list(table[1].str.strip())

print(f'Total entry: {len(names)}')

derepList = []
wgs, plm, psg, chm, cds = [[], [], [], [], []]
for name in names:
    isWgs, isPlm, isPsg, isChm, isCds = [False, False, False, False, False]
    if 'whole genome shotgun sequence' in name and 'plasmid' not in name:
        wgs.append(name)
        isWgs = True
    if 'plasmid' in name:
        if 'shotgun sequence' not in name:
            plm.append(name)
            isPlm = True
        else:
            psg.append(name)
            isPsg = True
    if ('complete genome' in name and 'plasmid' not in name) or \
        ('genome' in name and 'shotgun' not in name and 'complete genome' not in name and 'plasmid' not in name) or \
        ('chromosome' in name and 'whole genome shotgun sequence' not in name) or \
        ('chromosome' in name and 'complete sequence' in name) or \
            ('complete sequence' in name and 'plasmid' not in name):
        chm.append(name)
        isChm = True
    if 'complete CDS' in name:
        cds.append(name)
        isCds = True
    isKnown = sum([isWgs, isPlm, isPsg, isChm, isCds])
    if not isKnown:
        print(f'What is this: {name}?')
        break
    if isKnown > 1:
        print(name)
        print('Wgs' if isWgs else '')
        print('Plm' if isPlm else '')
        print('Psg' if isPsg else '')
        print('Chm' if isChm else '')
        print('Cds' if isCds else '')
        break

print(f'Total whole genome shotgun sequence {len(wgs)}')
print(f'Total plasmid {len(plm)}')
print(f'Total shotgun plasmid {len(psg)}')
print(f'Total chromosome {len(chm)}')
print(f'Total cds (single record) {len(cds)}')
if sum([len(wgs), len(plm), len(psg), len(chm), len(cds)]) != len(names):
    print('Unequal! Counted total:')
    print(sum([len(wgs), len(plm), len(psg), len(chm), len(cds)]))
    exit()
