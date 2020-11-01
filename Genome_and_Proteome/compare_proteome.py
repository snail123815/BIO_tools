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
Parser.add_argument("db", type=str, help="Target BLAST protein database path")
Parser.add_argument("q", type=str, help="Query protein file")
Parser.add_argument("-qname", type=str, help="name of input proteome", default="query")
Parser.add_argument("-tname", type=str, help="name of target proteome", default="target")
Parser.add_argument('-out', type=str, help="tsv format output file", default="")
Parser.add_argument("-eThresh", type=float, help="evalue threshold", default=0.001)
Parser.add_argument("-cThresh", type=int, help="coverage cutoff", default=60)

args = Parser.parse_args()

dbFile = args.db.strip()
protsFasta = args.q.strip()
qname = args.qname.strip()
tname = args.tname.strip()
tableOut = args.out.strip()
eThresh = args.eThresh
cThresh = args.cThresh

if tableOut == "":
    tableOut = os.path.splitext(protsFasta)[0] + f"_{tname}"
tableOut += "_corr"
tableOut += f"_e{eThresh}_cover{cThresh}.tsv"

dbPath, dbName = os.path.split(dbFile)
dbName = os.path.splitext(dbName)[0]

os.chdir(dbPath)

blastResTemp = NamedTemporaryFile()

blastpCmd = NcbiblastpCommandline(query=protsFasta, db=dbName, outfmt=5, out=blastResTemp.name, num_threads=os.cpu_count(),
                                  evalue=eThresh, max_hsps=1, max_target_seqs=1)
print("Running BLAST")
print(blastpCmd)

stdout, stderr = blastpCmd()

with open(tableOut, 'w') as outHandle:
    outHandle.write(f"{qname}\t{tname}\tCoverage\n")
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
