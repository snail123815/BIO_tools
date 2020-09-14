##################################################
# by Dù, Chāo (杜超)
# c.du@biology.leidenuniv.nl
# durand.dc@hotmail.com
# Institute of Biology Leiden, Leiden University
##################################################

'''
Usage:
$ python3 makeBLASTDB.py 
Sequence file(s) needs to be genbank or fasta format.
'''

from glob import glob
from subprocess import run
import os
import sys
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import tempfile
from time import sleep
import argparse


def steriliseName(name):
    """Remove all un-supported characters for BLAST program
    Do not pass a name with extension here
    """
    if name == None:
        return ""
    for t in '\\/:*?"<>|().':
        text = name.replace(t, ' ')
    return "_".join(text.split())


def getProteins(seqFile, nameTags=['locus_tag', 'protein_id', 'label']):
    """parse genebank file for protein sequences (all CDS)

    Arguments:
        seqFile {str} -- full path to genbank file

    Optional arguments:
        nameTags = ['locus_tag', 'protein_id']
            -- for a sequence of tags to extract protein names,
               will use the first tag found

    Returns:
        proteins -- list of protein SeqObjects
    """
    prots = []
    for seq in SeqIO.parse(seqFile, 'genbank'):
        for feat in seq.features:
            if feat.type == 'CDS':
                for nameTag in nameTags:
                    try:
                        protName = feat.qualifiers[nameTag][0]
                        break
                    except Exception as e:
                        if nameTag != nameTags[-1]:
                            continue  # not the last tag, try next
                        print(feat.qualifiers)  # for debug
                        raise e  # all tags not found
                try:
                    product = feat.qualifiers['product'][0]
                except:
                    product = ''
                try:
                    protId = feat.qualifiers['protein_id'][0]
                except:
                    protId = ''
                protSeq = feat.qualifiers['translation'][0]
                try:
                    prot = SeqRecord(
                        Seq(protSeq),
                        id=protName,
                        name=protId,
                        description=product
                    )
                except Exception as e:
                    print(protSeq)
                    print(protName)
                    print(product)
                    raise e
                prots.append(prot)
    return prots


def makedbGenbank(seqFile, outputPath):
    """make both nucl and protein blast databases from input genbank file
    Generate both nucleotides and protein fasta file
    Result database will appear in the same folder as the input file

    Arguments:
        seqFile {str} -- genbank full file path
                            fasta files of nucleotides and proteins will be stored 
                            beside the input file.

    Keyword Arguments:
        outputPath {str} -- destination to store database  

    Returns:
        outNuclDb, outProtDb -- full paths to database
    """
    # prepare file path and names
    seqPathName = os.path.splitext(seqFile)[0]
    seqName = os.path.split(seqPathName)[1]
    dbPath = outputPath
    if not os.path.isdir(dbPath):
        os.makedirs(dbPath)
        # makedirs support both simple folder and nested folder
        # mkdir only support simple folder
    seqName = steriliseName(seqName)
    nuclFile = f'{seqName}_nucl.fasta'
    protFile = f'{seqName}_prot.fasta'
    nuclPath = os.path.join(dbPath, nuclFile)
    protPath = os.path.join(dbPath, protFile)
    outNuclDb = os.path.join(dbPath, os.path.splitext(nuclFile)[0])
    outProtDb = os.path.join(dbPath, os.path.splitext(protFile)[0])

    # convert to nt fasta
    ntseqs = SeqIO.parse(seqFile, 'genbank')
    SeqIO.write(ntseqs, nuclPath, 'fasta')

    # parse CDS and write translation to prot fasta
    prots = getProteins(seqFile)
    SeqIO.write(prots, protPath, 'fasta')

    # make nt database
    runMakedbNucl = NcbimakeblastdbCommandline(
        dbtype='nucl',
        input_file=nuclPath,
        out=outNuclDb
    )
    stdout, stderr = runMakedbNucl()
    if stderr == '':
        print('\nSuccessfully made nucleotide database:')
        print(f'{outNuclDb}')
    else:
        print('Make nucleotide database failed')
        print('[stdout, stderr, seqPathName, seqName, dbPath, outNuclDb]')
        print('\n'.join([stdout, stderr, seqPathName,
                         seqName, dbPath, outNuclDb]))
        raise Exception()

    # make protein database
    runMakedbProt = NcbimakeblastdbCommandline(
        dbtype='prot',
        input_file=protPath,
        out=outProtDb
    )
    stdout, stderr = runMakedbProt()
    if stderr == '':
        print(
            f'\nSuccessfully made protein database\n{outProtDb}')
    else:
        print('Make protein database failed')
        print('[stdout, stderr, seqPathName, seqPath, seqName, dbPath, outProtDb]')
        print('\n'.join([stdout, stderr, seqPathName,
                         seqPath, seqName, dbPath, outProtDb]))
        raise Exception()

    return outNuclDb, outProtDb


def guessType(seqFile):
    """guess base on the first sequence of the fasta sequence

    Arguments:
        seqFile {str} -- full path to sequence file

    Returns:
        str -- 'nucl' or 'prot'
    """
    seqs = SeqIO.parse(seqFile, 'fasta')
    seqtype = 'nucl'
    firstSeq = str(next(seqs).seq).lower()  # only check the first sequence
    for n in 'atcgn*\n':
        # remove all linebreaks and possible nucleotide letters
        # ambiguous DNA is not supported in this simple function
        firstSeq = firstSeq.replace(n, '')
    if len(firstSeq) != 0:
        seqtype = 'prot'
    return seqtype


def makedbFasta(seqFile, outputPath):
    """make both nucl and protein database for blast from input file
    Generate both nucleotide and protein fasta file
    Result database will appear in the same folder as the input file

    Arguments:
        seqFile {str} -- input sequence (fasta)

    Return:
        dbtype, out -- dbtype = 'nucl' or 'prot
                       out = database path
    """
    # prepare path and file names
    seqPathName = os.path.splitext(seqFile)[0]
    seqName = os.path.split(seqPathName)[1]
    dbPath = outputPath
    if not os.path.isdir(dbPath):
        os.makedirs(dbPath)
        # makedirs support both simple folder and nested folder
        # mkdir only support simple folder
    seqName = steriliseName(seqName)
    out = os.path.join(dbPath, os.path.splitext(seqName)[0])
    dbtype = guessType(seqFile)
    # make database
    runMakedb = NcbimakeblastdbCommandline(
        dbtype=dbtype,
        input_file=seqFile,
        out=out
    )
    stdout, stderr = runMakedb()
    if stderr == '':
        print(
            f'\nSuccessfully made {dbtype} database\n{out}')
    else:
        print(f'Make {dbtype} database failed')
        print('[stdout, stderr, seqPathName, seqName, dbPath, out]')
        print('\n'.join([stdout, stderr, seqPathName,
                         seqName, dbPath, out]))
    return dbtype, out


def dbMade(seqFile, outputPath):
    """Check if a database from [seqFile] is made

    Arguments:
        seqFile {str} -- full path to a sequence (genbank/fasta) file

    Keyword Arguments:
        outputPath {str} -- relative path to the presumed output path (default: {'Blast_wrapper/GenomeProteomeDatabase'})

    Returns:
        bool -- True if the database is found, or the input do not have expected extensions.
                False otherwise. 
    """
    # prepare path and file names
    seqPathName, ext = os.path.splitext(seqFile)
    seqPath, seqName = os.path.split(seqPathName)
    dbPath = os.path.join(seqPath, outputPath)
    seqName = steriliseName(seqName)
    # store possible db files in list
    if ext in ['.fa', '.fasta', '.faa']:
        dbs = [os.path.join(dbPath, f) for f in [f'{seqName}.psq',
                                                 f'{seqName}.nsq']]
        # Either protein or nucleotide database does not matter
        # as the info is in seqName
    elif ext in ['.gb', '.genbank', '.gbk']:
        dbs = [os.path.join(dbPath, f) for f in [f'{seqName}_prot.psq',
                                                 f'{seqName}_nucl.nsq']]
    else:  # file extension do not match
        return True
    # check if any possible db file exists
    for db in dbs:
        if os.path.isfile(db):
            print(f'\nThis database is made, skip: {db}')
            return True
    return False


def makedb(f, outputPath):
    """Make blast database(s) for input file f
    Try to make both protein and nucleotide database for genbank file,
    determin database type automatically from fasta file.

    Arguments:
        f {str} -- full path to the sequence file (genbank/fasta)

    Keyword Arguments:
        outputPath {str} -- relative path to the database (relative to the input sequence)
                               (default: {'Blast_wrapper/GenomeProteomeDatabase'})

    Returns:
        outNuclDb, outProtDb -- paths to the generated databases
    """
    ext = os.path.splitext(f)[1]
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)
    outNuclDb, outProtDb = ('', '')
    if dbMade(f, outputPath):
        pass
    else:
        if ext in ['.gb', '.genbank', '.gbk']:
            outNuclDb, outProtDb = makedbGenbank(f, outputPath)
        elif ext in ['.fa', '.fasta', '.faa']:
            dbtype, out = makedbFasta(f, outputPath)
            if dbtype == 'nucl':
                outNuclDb = out
            else:
                outProtDb = out
    return outNuclDb, outProtDb


def combinDbs(dbList, dbtype, dbName):
    """combin databases in dbList to a database called dbName

    Arguments:
        dbList {list} -- list of paths to database
        dbtype {str} -- 'nucl' or 'prot'
        dbName {str} -- output database name
    """
    # prepare a file containing all database paths
    # I find it does not work if you pass database names to the command line
    # directly by "' '.join(list)" string
    listFile = tempfile.NamedTemporaryFile('w')
    for db in dbList:
        listFile.write(db)
        listFile.write('\n')
    listFile.seek(0)  # this will actually write the file to disk
    # generate an output path using the path of the first database
    outputDb = os.path.join(os.path.split(dbList[0])[0], f'{dbName}_{dbtype}')
    combinDbArgs = [
        'blastdb_aliastool',
        '-dblist_file', listFile.name,
        '-dbtype', dbtype,
        '-out', outputDb,
        '-title', f'{dbName}_{dbtype}'
    ]
    combinDb = run(combinDbArgs, capture_output=True)
    if combinDb.returncode != 0:
        print()
        print(' '.join(combinDb.args))
        print(combinDb.stdout.decode())
        print(combinDb.stderr.decode())
        raise Exception('Combine db failed')
    print(f'\n{dbtype} db combined:\n{outputDb}')


def prepareArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input', help="path to sequence file or folder with sequence files", nargs='+')
    parser.add_argument('-t', '--targetPath',
                        help='target path to store generated database')
    parser.add_argument('-n', '--databaseName', help='''
    If multiple database was made, the program try to combine them, 
    this parameter specifies the name of combined database, will be sterialized.
    ''')
    args = parser.parse_args()
    ipaths = args.input
    knownExts = ['.fasta', '.fa', '.faa', '.gbk', '.gbff', '.gb']
    filePaths = []
    for p in ipaths:
        if os.path.isdir(p):
            for ext in knownExts:
                filePaths += glob(os.path.join(p, f'*{ext}'))
        else:
            if os.path.splitext(p)[1] not in knownExts:
                sys.exit(f'Extension not known ({" ".join(knownExts)})')
            filePaths.append(p)
    if args.targetPath == None:
        dbpath = os.path.split((os.path.splitext(ipaths[0])[0]))[
            0]  # path where dir or file located
    else:
        dbpath = args.t.strip()
    dbpath = os.path.join(dbpath, 'blastdb')
    if args.databaseName != None:
        if len(filePaths) == 1:
            sys.exit('only one sequence, no need to combine database')
        combinedDbname = steriliseName(args.n.strip())
    else:
        combinedDbname = os.path.split((os.path.splitext(ipaths[0])[0]))[1]
    return filePaths, dbpath, combinedDbname


def main(filePaths, dbpath, combinedDbname):
    # initialize lists for output databases
    nuclDbs = []
    protDbs = []
    # make database
    for f in filePaths:
        outNuclDb, outProtDb = makedb(f, outputPath=dbpath)
        if outNuclDb != '':
            nuclDbs.append(outNuclDb)
        if outProtDb != '':
            protDbs.append(outProtDb)
    # combine database if possible
    for dbList, dbtype in zip([nuclDbs, protDbs], ['nucl', 'prot']):
        if len(dbList) <= 1:  # no need to combine
            continue
        combinDbs(dbList, dbtype, combinedDbname)


if __name__ == '__main__':
    filePaths, dbpath, combinedDbname = prepareArguments()
    main(filePaths, dbpath, combinedDbname)
