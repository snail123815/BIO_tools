from subprocess import run
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os
import sys
import platform
from tempfile import NamedTemporaryFile

# For CRISPR spacer design: generate seqs for BLAST checking the specificity.
# Add blast feature later
cwd = os.path.dirname(os.path.realpath(__file__))

# Define path


def defAppPath():
    if platform.system() == 'Darwin':
        blastn_program_path = 'blastn'
        fileReader = ['open', '-a', 'TextEdit']
    elif platform.system() == 'Linux':
        blastn_program_path = 'blastn'
        fileReader = False
    else:  # 'Windows'
        blastn_program_path = r"D:\Program Files\NCBI\blast-2.6.0+\bin\blastn.exe"
        fileReader = ['notepad']
    return blastn_program_path, fileReader


def generate_4seq(seq):
    A_seq = SeqRecord(seq + Seq('AGG', IUPAC.unambiguous_dna),
                      'A', description='')
    T_seq = SeqRecord(seq + Seq('TGG', IUPAC.unambiguous_dna),
                      'T', description='')
    C_seq = SeqRecord(seq + Seq('CGG', IUPAC.unambiguous_dna),
                      'C', description='')
    G_seq = SeqRecord(seq + Seq('GGG', IUPAC.unambiguous_dna),
                      'G', description='')
    all_sps = []
    all_sps.append(A_seq)
    all_sps.append(T_seq)
    all_sps.append(C_seq)
    all_sps.append(G_seq)
    return all_sps


def getSeqInput():
    my_seq = Seq(input('23 bp, spacer + PAM sequence:\n'),
                 IUPAC.unambiguous_dna)
    if len(my_seq) == 20:
        print(f'\nSpacer sequence:\n{my_seq}\nOrigional PAM unknown')
        rc_sp = my_seq.reverse_complement()
        my_seq_core = my_seq[8:]
        print(f'Core sequence: {my_seq_core}+PAM\n')
        print(f'Reverse complemented spacer:\n{rc_sp}\n')
        all_sps = generate_4seq(my_seq)
    elif len(my_seq) == 23:
        print(
            f'\nSpacer sequence:\n{my_seq[:-3]}\nOrigional PAM = {my_seq[-3:]}')
        rc_sp = my_seq[:-3].reverse_complement()
        my_seq_core = my_seq[8:-3]
        print(f'Core sequence: {my_seq_core}+PAM\n')
        print(f'Reverse complemented spacer:\n{rc_sp}\n')
        all_sps = generate_4seq(my_seq[:-3])
    else:
        raise Exception('Wrong input!')
    return my_seq, rc_sp, all_sps


def printDesign(my_seq, rc_sp, all_sps):
    print('\nFor pCRISPomyces-2 one spacer design:')
    print(f'Forward spacer: acgc{my_seq[:-3].upper()}')
    print(f'Reverse spacer: aaac{rc_sp.upper()}\n')

    print('\nFor pCRISPR-Cas9 and pCRISPR-dCas9 design:')
    print(f'catgccatgg{my_seq[:-3].upper()}gttttagagctagaaatagc')


def setBlastDBEnv(dbPath=None, origionalBLASTDBenv=None):
    if dbPath == None:
        print('\nRestoring environment varialbe... ', end='')
        if origionalBLASTDBenv:
            os.environ['BLASTDB'] = origionalBLASTDBenv
        else:
            os.environ.pop('BLASTDB')
        print('Restored!')
        return

    print('\nSetting up BLASTDB environment variable...')
    if "BLASTDB" in os.environ:
        origionalBLASTDBenv = os.environ('BLASTDB')
    else:
        origionalBLASTDBenv = False
    try:
        os.environ['BLASTDB'] = dbPath
        print(dbPath)
        print('Success, environment variable BLASTDB will be restored in the end.\n')
    except:
        print('BLASTDB environment variable set up failed.')
        input('<Enter> to exit')
        exit()
    return origionalBLASTDBenv


def doBlast(all_sps, nuclDb, blastn_program_path='blastn', fileReader=False):
    tempInput = NamedTemporaryFile()
    tempOutput = NamedTemporaryFile()
    dbPath, db = os.path.split(nuclDb)
    if os.path.splitext(db)[1] in [
        '.nal', '.pal', '.ndb', '.nhr', '.nin', '.not', '.nsq', 
        '.ntf', '.nto', '.phr', '.pin', '.pot', '.psq', '.ptf', '.pto'
        ]:
        db = os.path.splitext(db)[0]
    # write spacers to tempfile
    SeqIO.write(all_sps, tempInput.name, 'fasta')

    num_threads = os.cpu_count()
    run_blastn = NcbiblastnCommandline(
        blastn_program_path,
        query=tempInput.name,
        out=tempOutput.name,
        db=db,
        remote=False,
        evalue=0.1,
        num_threads=num_threads,
        word_size=6
    )
    # set BLASTDB environment
    origionalBLASTDBenv = setBlastDBEnv(dbPath=dbPath)
    try:
        stdout, stderr = run_blastn()
        if stderr == '':
            print(f'\nSuccess!!\n{tempOutput.name}\n')
        else:
            print(stdout.decode())
            print(stderr.decode())
    except Exception as e:
        print(e)
    finally:
        setBlastDBEnv(origionalBLASTDBenv=origionalBLASTDBenv)
    tempOutput.seek(0)
    if fileReader != False:
        try:
            run([*fileReader, tempOutput.name])
        except Exception as e:
            print(e)
            print(tempOutput.readlines())


def main():
    blastn_program_path, fileReader = defAppPath()
    nuclDb = sys.argv[1]

    while True:
        my_seq, rc_sp, all_sps = getSeqInput()
        printDesign(my_seq, rc_sp, all_sps)
        doBlast(all_sps, nuclDb, blastn_program_path, fileReader=fileReader)


if __name__ == '__main__':
    main()
