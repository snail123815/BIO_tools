#######################################################
# by Chao DU
# c.du@biology.leidenuniv.nl
#######################################################

import os
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def simplifySeqRec(rec):
    newRec = SeqRecord(
        rec.seq,
        id=rec.id,
        name=rec.name,
        description=rec.description,
        dbxrefs=rec.dbxrefs,
        annotations=rec.annotations,
        letter_annotations=rec.letter_annotations
    )
    for feat in rec.features:
        if feat.type not in ['CDS', 'gene', 'source']:
            continue
        newRec.features.append(feat)
    return newRec

def simplifyGBK(file, out):
    newRecs = []
    records = SeqIO.parse(file, 'genbank')
    for rec in records:
        newRecs.append(simplifySeqRec(rec))
    SeqIO.write(newRecs, out, 'genbank')
    return


def main():
    aparser = argparse.ArgumentParser(description="Keep only CDS, gene, and source")
    aparser.add_argument('path', help="path to gbk or dir of gbks")

    args = aparser.parse_args()
    inputPath = args.path
    files = []
    outFiles = []
    if os.path.isdir(inputPath):
        outdir = f'{inputPath}_simple'
        while True:
            try:
                os.mkdir(outdir)
                break
            except FileExistsError:
                outdir = f'{outdir}_x'
        for _, _, fs in os.walk(inputPath):
            files = [os.path.join(inputPath, f) for f in fs]
            for f in fs:
                fn, ext = os.path.splitext(f)
                outFiles.append(os.path.join(outdir, f'{fn}_simple{ext}'))
            break
    else:
        fn, ext = os.path.splitext(inputPath)
        files.append(inputPath)
        outFiles.append(f'{fn}_simple{ext}')
    for f, newf in zip(files, outFiles):
        simplifyGBK(f, newf)
    pass


if __name__ == "__main__":
    main()
