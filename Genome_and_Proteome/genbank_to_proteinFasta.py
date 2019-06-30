from Bio import SeqIO

genbankSeq = SeqIO.parse(
    '/Users/durand.dc/Documents/works/Resources/Strains and Plasmids/MBT84.gbk', 'genbank')

proteins = {}

for seq in genbankSeq:
    for feat in seq.features:
        if feat.type == 'CDS':
            gene = feat.qualifiers['gene'][0].capitalize(
            ) if 'gene' in feat.qualifiers else None
            locusTag = feat.qualifiers['locus_tag'][0] if 'locus_tag' in feat.qualifiers else None
            proteinName = '_'.join(filter(None, (locusTag, gene))).replace(' ', '_')
            proteins[proteinName] = feat.qualifiers['translation'][0]

# for feat in genbankSeq.features:
#     hasTranslationQualifier = 'translation' in feat.qualifiers
#     hasTranslationNote = False
#     if not hasTranslationQualifier:
#         if 'note' in feat.qualifiers:
#             for note in feat.qualifiers['note']:
#                 if note.startswith('/translation'):
#                     hasTranslationNote = True
#                     translation = note.split('=')[-1]
#     else:
#         translation = feat.qualifiers['translation'][0]
#
#     if hasTranslationQualifier or hasTranslationNote:
#         prefix = 'pET28a'
#         proteinId = feat.qualifiers['protein_id'][0] if 'protein_id' in feat.qualifiers else None
#         gene = feat.qualifiers['gene'][0].capitalize() if 'gene' in feat.qualifiers else None
#         label = feat.qualifiers['label'][0] if 'label' in feat.qualifiers and not gene else None
#         proteinName = '_'.join(filter(None, (prefix,proteinId,gene,label))).replace(' ', '_')
#         proteins[proteinName] = translation

with open('/Users/durand.dc/Documents/works/Resources/Strains and Plasmids/MBT84.proteome.fasta', 'w') as output:
    for protein in proteins:
        output.write(f'>{protein}\n')
        output.write(proteins[protein])
        output.write('\n')
