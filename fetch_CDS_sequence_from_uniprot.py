# Fetch CDS DNA sequence from only Uniprot accession

# First, it is not accurate!! Because one uniprot protein often comes from multiple CDS sequences, some times they are identical, but some times they are not!

# Second, the reference records are fetched swissprot records


# re is for get number part from text. This is because sometimes the returned position is like '>4232312'.
import re
from math import ceil
from Bio import Entrez
from Bio import ExPASy
from Bio import SwissProt

for accession in accessionBin:
    records = []
    handle = ExPASy.get_sprot_raw(accession)
    try:
        record = SwissProt.read(handle)
    except ValueException:
        print(f"WARNING: Accession {accession} not found")
        continue
    uniprotRes[bin].append(record)


# Refined version

# Prepare Entrez accession list from uniprot records


Entrez.email = 'c.du@biology.leidenuniv.nl'

# The origional record name stored in a dict
badRequests = []
print(f'Total {len(uniprotRec)}')
emblAccs = {}
for record in uniprotRec:
    if record.cross_references[0][0] != 'EMBL':
        raise  # it should alwasy be 'EMBL'
    emblAcc = record.cross_references[0][2]
    emblAcc = emblAcc.split('.')[0]
    emblAccs[emblAcc] = record.entry_name
# get protein in genbank format
allProteinGb = []
accsList = list(emblAccs.keys())

# Fetch protein data in groups
# The nucleotide ID and the CDS position for the protein can be found in the gb file fetched


# Fetch data in groups will make the process faster. Upper limit suggested by NCBI is 200.
# Here I do 150

for i in range(ceil(len(emblAccs) / 150)):
    print(i, end='|')
    if i < len(emblAccs) // 150:
        accs = accsList[i * 150:(i + 1) * 150]
    else:
        accs = accsList[i * 150:]
    handle = Entrez.efetch(db="protein", id=accs, rettype='gb', retmode="text")
    records = SeqIO.parse(handle, 'genbank')
    for rec in records:
        allProteinGb.append(rec)
    handle.close()

print(len(allProteinGb))

# Parse protein record
# Get the nucleotide ID and CDS position, store in a dict


cdsDb = {}  # to store the ID and position
print(f'Total {len(allProteinGb)}')
count = 0
for rec in allProteinGb:
    count += 1
    print(count, end='|')
    for feat in rec.features:
        if feat.type != 'CDS':
            continue
        if 'coded_by' not in feat.qualifiers:
            print(feat.qualifiers)
            break
        loctext = feat.qualifiers.get('coded_by')[0]
        if loctext.startswith('complement'):
            strand = 2
            loctext = loctext.split('(')[1][:-1]
        else:
            strand = 1

        nuclId, locations = loctext.split(':')
        nuclId = nuclId.split('.')[0]
        start, stop = locations.split('..')
        start = re.findall(r'\d+', start)[0]
        stop = re.findall(r'\d+', stop)[0]
        cdsDb[nuclId] = [emblAccs[rec.name], strand, start, stop]
        break


# Fetch DNA sequence
# The sequence has to be fetched one by one because I don't know the way to pass locations in general

print(f'Total {len(cdsDb)}')
ignoredSeqs = []
seqs = []
count = 0
for nuclId in cdsDb:
    count += 1
    print(f'{count}', end='|')
    entry_name, strand, start, stop = cdsDb[nuclId]
    handle = Entrez.efetch(db='nuccore', id=nuclId, rettype='fasta',
                           retmode='text', strand=strand, seq_start=start, seq_stop=stop)
    seq = SeqIO.read(handle, 'fasta')
    seq.name = f"{emblAcc}_CDS"
    seq.id = entry_name
    seq.description = f'{emblAcc}_CDS'
    if len(seq) > 300:
        print(f'{count}ignored', end='|')
        ignoredSeqs.append(seq)
    else:
        seqs.append(seq)
keptCdsFile = 'CDS_DNAseqs.txt'
ignoredSeqsFile = 'longer_DNAseqs.txt'
print(f'Writing {len(seqs)} seqs to {keptCdsFile}...')
SeqIO.write(seqs, keptCdsFile, 'fasta')
print(f'Writing {len(ignoredSeqs)} seqs to {keptCdsFile}...')
SeqIO.write(ignoredSeqs, ignoredSeqsFile, 'fasta')
