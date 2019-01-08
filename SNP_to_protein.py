import os
from Bio import SeqIO
from Bio.Seq import reverse_complement, Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from subprocess import Popen, PIPE, STDOUT

# read the proteome file to change
baseProteome_path = '/Users/durand.dc/Documents/works/Resources/Resource_M145/M145+plasmids_proteome_combined_strepDBandgenbank_20181011.fasta'
baseProteome = {}
for seq in SeqIO.parse(baseProteome_path, 'fasta'):
    id = seq.id.split(' ')[0]
    baseProteome[id] = seq
# print(baseProteome['SCO1839'])

# read genome file
genomePath = '/Users/durand.dc/Documents/works/Resources/Resource_M145/M145PlusStredbNotes.gb'
genome = SeqIO.read(genomePath, 'genbank')
locationDict = {}
print('Parsing genome...')
for feat in genome.features:
    if feat.type == 'CDS':
        if 'gene' in feat.qualifiers:
            gene = feat.qualifiers['gene'][0].upper()
        else:
            gene = ''
        if 'locus_tag' in feat.qualifiers:
            locusTag = feat.qualifiers['locus_tag'][0].upper()
        else:
            locusTag = ''
        # only SCO is considered:
        if gene.startswith('SCO'):
            id = gene
        elif locusTag.startswith('SCO'):
            id = locusTag
        else:
            id = '_'.join(filter(None, (locusTag, gene))).replace(' ', '_')
        if len(id) != 7:
            print(id)
        locationDict[id] = feat.location
print('Above are non-uniform IDs.\n')
# print(list(locationDict['SCO1839']))

# read the change table
snpTable_path = '/Users/durand.dc/Documents/works/Project_SYSTERACT/SNPs_genome_variants/SYSTERACT_SNP_coding.xlsx'
snpTable = pd.read_excel(snpTable_path, index_col='Position')


def printChanges(a, b, clustalFilePath=''):
    proc = Popen(['clustalo', '-i', '-', '--outfmt', 'clu', '--residuenumber'],
                 stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    input = f'>{a.id}\n{a.seq}\n>{b.id}\n{b.seq}'
    stdout = proc.communicate(input=bytes(input, 'utf-8'))[0]
    if clustalFilePath == '':
        print(stdout.decode('utf-8'))
    else:
        with open(clustalFilePath, 'a') as handle:
            handle.write(stdout.decode('utf-8'))
            handle.write(f'\n{"-"*100}\n')
# printChanges


def mutProt(protId, positions, oris, vars, modes, baseProteome=baseProteome, locationDict=locationDict, genome=genome, clustalFilePath=''):
    """mode = ['change', 'dup', 'del', 'ins']"""
    if protId not in baseProteome.keys():
        print(f'{protId} not in base proteome')
        return
    if protId not in locationDict.keys():
        print(f'{protId} not in genome')
        return
    baseProtein = baseProteome[protId]
    start = locationDict[protId].start
    end = locationDict[protId].end
    if locationDict[protId].strand == -1:
        oris = [reverse_complement(ori) for ori in oris]
        vars = [reverse_complement(var) for var in vars]
    modifiedSeq = genome.seq.tomutable()
    for i, pos in enumerate(positions):
        pos = int(pos)
        mode = modes[i]
        ori = oris[i]
        var = vars[i]
        if mode in ['change', 'del']:
            oriFromPos = str(genome[pos - 1:pos - 1 + len(ori)].seq)
            if oriFromPos != ori:
                print(protId)
                print([pos - 1, pos - 1 + len(ori)])
                print(f'ori From Pos:{oriFromPos}')
                print(f'ori Form Input:{ori}')
                raise
        if mode == 'change':
            modifiedSeq[pos - 1] = var
        if mode in ['dup', 'ins']:
            modifiedSeq = modifiedSeq[:pos] + var + modifiedSeq[pos:]
            end += len(var)
        if mode == 'del':
            modifiedSeq = modifiedSeq[:pos] + modifiedSeq[pos + len(ori):]
            end -= len(ori)
    modifiedSeq = modifiedSeq.toseq()

    if locationDict[protId].strand == 1:
        newSeqPart = modifiedSeq[start:]
    else:  # reverse complement
        newSeqPart = modifiedSeq[:end].reverse_complement()
    endTri = len(newSeqPart) // 3 * 3
    newSeqPart = newSeqPart[:endTri]
    newProtSeq = newSeqPart.translate(to_stop=True).tomutable()
    newProtSeq[0] = 'M'
    id = f'{baseProtein.id}_{strainMutCol}'
    description = baseProtein.description.replace(f'{baseProtein.id} ', "")
    newProt = SeqRecord(newProtSeq, id=id,
                        description=f'{description} variant {strainMutCol}')
    printChanges(baseProtein, newProt, clustalFilePath=clustalFilePath)
    return newProt
# mutProt


def getMut(protein, clustalFilePath=''):
    snpTable_prot = snpTable[snpTable.loc[:, strainMutCol] == protein]
    positions = []
    modes = []
    oris = []
    vars = []
    for i, c in enumerate(list(snpTable_prot.loc[:, 'HGVS.c'])):
        if '>' in c:
            position = snpTable_prot.index[i]
            mode = 'change'
            ori, var = c.split('>')
            ori = ori[-1]
        elif 'dup' in c:
            position = snpTable_prot.index[i]
            mode = 'dup'
            ori = ''
            var = c.split('dup')[1]
        elif 'ins' in c:
            position = snpTable_prot.index[i]
            mode = 'ins'
            ori = ''
            var = c.split('ins')[1]
        elif 'del' in c:
            position = snpTable_prot.index[i]
            mode = 'del'
            ori = c.split('del')[1]
            var = ''
        else:
            raise
        positions.append(position)
        modes.append(mode)
        oris.append(ori)
        vars.append(var)
    return mutProt(protein, positions, oris, vars, modes, clustalFilePath=clustalFilePath)
# getMut


for strainMutCol in ['mut_M1152', 'mut_M145&MatAB', 'mut_common']:
    clustalFilePath = f'/Users/durand.dc/Documents/works/Project_SYSTERACT/SNPs_genome_variants/SYSTERACT_SNP_coding_{strainMutCol}.txt'
    if os.path.isfile(clustalFilePath):
        os.remove(clustalFilePath)
    allProteins = snpTable.loc[:, strainMutCol].drop_duplicates().dropna()
    newProteome = baseProteome.copy()
    for protein in allProteins:
        print(protein)
        newProt = getMut(protein, clustalFilePath=clustalFilePath)
        del newProteome[protein]
        if len(newProt.seq) < 20:
            continue
        else:
            newProteome[protein] = newProt

    SeqIO.write(newProteome.values(),
                f'/Users/durand.dc/Documents/works/Project_SYSTERACT/SNPs_genome_variants/M145+plasmids_proteome_combined_strepDBandgenbank_RNASeq_var_{strainMutCol}.fasta', 'fasta')
