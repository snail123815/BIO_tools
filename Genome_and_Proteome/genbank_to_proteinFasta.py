from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def getProteins(seq):
    """Extract proteins from a SeqRecord. If your input file have multiple contigs, do a loop"""
    proteins = []
    for feat in seq.features:
        if feat.type == 'CDS':
            gene = feat.qualifiers['gene'][0] if 'gene' in feat.qualifiers else None
            # make upper case for this is a protein
            if gene != None:
                gene = gene[0].upper() + gene[1:]
            locusTag = feat.qualifiers['locus_tag'][0] if 'locus_tag' in feat.qualifiers else None
            proteinName = ', '.join(filter(None, (gene, locusTag))).replace(' ', '_')
            prot = SeqRecord(Seq(feat.qualifiers['translation'][0]),
                             id=locusTag, description=proteinName)
            proteins.append(prot)
    return proteins


if __name__ == "__main__":
    import os
    import argparse
    myParser = argparse.ArgumentParser(description='Extract CDS from genbank file')
    myParser.add_argument('path_gbk', help='path to your genbank file')
    args = myParser.parse_args()
    input = args.path_gbk
    output = f'{os.path.splitext(input)[0]}_proteins.fa'
    seqIn = SeqIO.parse(input, 'genbank')

    proteins = []
    for seq in seqIn:
        proteins.extend(getProteins(seq))
    print(f"Total CDS {len(proteins)}")

    SeqIO.write(proteins, output, 'fasta')
    print(f'Done, check file\n{output}')
