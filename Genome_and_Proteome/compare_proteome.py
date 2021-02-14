about = """
Build correlation between two genomes, find homologue for each protein.
BLAST every protein in query file against target database,
give the best hit that score (e and coverage) higher than required threshold.

By DU Chao
c.du@biology.leidenuniv.nl
durand.dc@hotmail.com
"""

from os import ctermid
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import argparse
import os
from tempfile import NamedTemporaryFile


Parser = argparse.ArgumentParser(description=about, formatter_class=argparse.RawDescriptionHelpFormatter)
Parser.add_argument("targetDatabase", type=str, help="Target protein BLAST database path")
Parser.add_argument("queryProteinFasta", type=str, help="Query protein file in fasta format")
Parser.add_argument("-qname", type=str, help="name of input proteome", default="query")
Parser.add_argument("-tDbName", type=str, help="name of target proteome", default="name of target database")
Parser.add_argument('-out', type=str, help="tsv format output file, default to [queryFile]_[tDbName]_corr_e[evalue]_cover[coverage].tsv", default="")
Parser.add_argument("-eThresh", type=float, help="evalue threshold", default=0.001)
Parser.add_argument("-cThresh", type=int, help="coverage cutoff", default=60)

args = Parser.parse_args()

dbFile = args.targetDatabase.strip()
protsFasta = args.queryProteinFasta.strip()
qname = args.qname.strip()
tDbName = args.tDbName.strip()
tableOut = args.out.strip()
eThresh = args.eThresh
cThresh = args.cThresh

for n in [dbFile, protsFasta, qname, tDbName, tableOut]:
    fn = os.path.splitext(n)[0]
    for s in ['.', ' ']:
        assert s not in fn, f'"{s}" not allowed in the file names. Check: {n}.'


if tableOut == "":
    tableOut = os.path.splitext(protsFasta)[0] + f"_{tDbName}"
tableOut += "_corr"
tableOut += f"_e{eThresh}_cover{cThresh}.tsv"

dbPath, dbName = os.path.split(dbFile)
dbName = os.path.splitext(dbName)[0]

if dbPath != '': # current dir
    os.chdir(dbPath)

blastResTemp = NamedTemporaryFile()

blastpCmd = NcbiblastpCommandline(query=protsFasta, db=dbName, outfmt=5, out=blastResTemp.name, num_threads=os.cpu_count(),
                                  evalue=eThresh, max_hsps=1, max_target_seqs=1)
print("Running BLAST")
print(blastpCmd)

stdout, stderr = blastpCmd()

print(f'Writing result in table {tableOut}')
with open(tableOut, 'w') as outHandle:
    outHandle.write(f"{qname}\t{tDbName}\tCoverage\n")
    with open(blastResTemp.name, 'r') as resHandle:
        records = NCBIXML.parse(resHandle)
        for record in records:
            queryName = record.query
            target = "NotFound"
            coverage = 0
            if len(record.alignments) != 0:
                firstAlignment = record.alignments[0]
                hsp = firstAlignment.hsps[0]
                hspCoverage = hsp.align_length/record.query_length
                if hspCoverage >= cThresh/100:
                    target = firstAlignment.hit_def
                    coverage = hspCoverage
            outHandle.write(f"{queryName}\t{target}\t{coverage:.0%}\n")
