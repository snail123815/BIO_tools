
import argparse
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def getProteins(s):
    ps = []
    for feat in s.features:
        if feat.type == "CDS":
            pid = feat.qualifiers['locus_tag'][0]
            try:
                pgid = feat.qualifiers['gene'][0]
                pgid = pgid[0].upper()+pgid[1:]
            except KeyError:
                pgid = ''
            try:
                pseq = Seq(feat.qualifiers['translation'][0])
            except KeyError:
                pgene = s.seq[feat.location.start:feat.location.end]
                pgene = (pgene if feat.location.strand == 1 else pgene.reverse_complement())
                pseq = pgene.translate()
            pfun = feat.qualifiers['product'][0]
            p = SeqRecord(pseq, id=pid, name=pgid, description=pfun)
            ps.append(p)
    return ps


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('file', help='genbank file')

    args = argparser.parse_args()
    filePath = os.path.realpath(args.file)
    outputPath = filePath+'.prot.fa'
    inputSeqs = list(SeqIO.parse(filePath, 'genbank'))

    ps = []
    for s in inputSeqs:
        ps.extend(getProteins(s))

    n = SeqIO.write(ps, outputPath, 'fasta')
    print(f'Successfully wrote {n} proteins')

if __name__ == "__main__":
    main()