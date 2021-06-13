#########################################
# By: Du, Chao
# c.du@biology.leidenuniv.nl
# MBT, IBL, Leiden Univ.
# April 2021
##########################################

blastpFetchCDSAlign_description='''
    This script does following:
    1. Do protein BLAST for input file on NCBI server
    2. Find identical proteins for all results and Fetch coding sequence
       Only get one identical protein from one organism.
    3. Align these sequences use `clustalo`
    '''
blastpFetchCDSAlign_epilog='Prerequisites: blastp and clustalo.'

import argparse
import os
import time
from subprocess import run

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import Entrez


def main():
    parser = argparse.ArgumentParser(description=blastpFetchCDSAlign_description,
                                     epilog=blastpFetchCDSAlign_epilog)
    parser.add_argument('email', help='Email address for Entrez database query')
    parser.add_argument('api_key', help="NCBI's API Key for Entrez database query, you can get it from 'Account settings' of your NCBI account")
    parser.add_argument('fasta', help='Protein fasta file')
    parser.add_argument('--ident', type=float, help='Min identity percentage (hsp.identities/len(hsp.query)) to use',
                        default=0.4)

    args = parser.parse_args()


    fastaFile = args.fasta
    Entrez.api_key = args.api_key
    Entrez.email = args.email
    basePath, name = os.path.split(os.path.realpath(fastaFile))
    name = os.path.splitext(name)[0]
    identT = args.ident

    print()
    print('JOB DISCRIPTION')
    print(parser.description)
    print('Entrez ID - email:     '+Entrez.email)
    print('Entrez ID - API key:   '+Entrez.api_key)
    print('Query file:            '+fastaFile)
    print('Threshold of identity: '+f'{identT:.0%}')

    # Do BLAST
    blastpResXML = os.path.join(basePath, f'{name}.blastp.xml')
    startTime = time.time()
    if not os.path.isfile(blastpResXML):
        print()
        print('-'*80,'\nDoing blastp.')
        cline = NcbiblastpCommandline(query=fastaFile, db='nr', remote=True, out=blastpResXML, outfmt=5)
        print(cline)
        cline()
        print(f"Blast time: {time.time()-startTime:.2f}")
        print('Done.')
    else:
        print()
        print('-'*80+f'\nUse existing blastp result {blastpResXML}')
    with open(blastpResXML, 'rb') as xmlhandle:
        records = list(NCBIXML.parse(xmlhandle))

    # Parse BLAST result and fetch CDS from database
    print()
    print('-'*80,'\nSearch for identical proteins and fetch CDS:')
    resultFils = []
    for i, rec in enumerate(records):
        fetchingStartTime = time.time()
        print(f'Record {i+1}/{len(records)}')
        id = rec.query.replace(' ', '_').replace('/', '_')
        cdsFasta = os.path.join(basePath, f'{id}_BLASTPres_corrCDSs.fasta')
        resultFils.append(cdsFasta)
        print('CDS file: ', cdsFasta)
        if os.path.isfile(cdsFasta):
            print('Use existing cds file.')
            continue
        resultAccs = getBlastpProteins(rec, identT)
        cdsList = []  # SeqRecords
        nresults = len(resultAccs)
        for x, acc in enumerate(resultAccs):
            print(f'\r{x+1}/{nresults}', end=('' if x+1 < nresults else '\n'))
            cdsList.extend(getCdsFromProtein(acc))
        SeqIO.write(cdsList, cdsFasta, 'fasta')
        print(f'Fetching time: {time.time() - fetchingStartTime:.2f}')
    print('Done.')

    # Run Alignment
    print()
    print('-'*80,'\nRunning alignments:')
    alignmentStartTime = time.time()
    for i, f in enumerate(resultFils):
        print(f'Record {i+1}/{len(resultFils)}')
        alignOut = f'{os.path.splitext(f)[0]}.alignment.clu'
        args = [
            'clustalo',
            '-i', f,
            '-o', alignOut,
            '--outfmt', 'clu',
            '--threads', f'{os.cpu_count()}',
            '--force'
        ]
        print(alignOut)
        if not os.path.isfile(alignOut):
            print(' '.join(args))
            run(args)
            print(f'Time elapsed since alignment start: {time.time() - alignmentStartTime:.2f}')
    print()
    print('-'*80,f'\nTotal time used: {time.time() - startTime:.2f}')
    print('Done')


def getBlastpProteins(rec, identT):
    resultAccs = []
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            acc = alignment.accession
            identprec = hsp.identities/len(hsp.query)
            if identprec >= identT:
                resultAccs.append(acc)
    resultAccs = list(set(resultAccs))
    resultAccs.sort()
    return resultAccs


def getCdsFromProtein(acc):
    # Process every accession one by one rather than
    # request multiple accessions at once. It seems Entrez.read()
    # can only read the last result, and the return value of
    # current Entrez is not readable by Entrez.parse().
    # Besides, when fetching coding sequences, it is only
    # feasible to do it one by one with a start and stop position.
    organisms = []
    resultCdsList = []
    # First query the identical proteins database.
    with Entrez.efetch(db='protein', rettype='ipg', retmode='xml', id=acc) as handle:
        ipgRec = Entrez.read(handle)
    ipList = ipgRec['IPGReport']['ProteinList']
    for protRec in ipList:
        # For each query, we need to see if it corresponds to different organism.
        # There is chance that identical proteins are encoded differently.
        recAcc = protRec.attributes['accver']
        cds = protRec['CDSList'][0]
        nucAcc = cds.attributes['accver']
        start = cds.attributes['start']
        stop = cds.attributes['stop']
        strand = cds.attributes['strand']
        org = cds.attributes['org']
        if org in organisms:
            continue
            # Use only one protein from single organism
        organisms.append(org)
        with Entrez.efetch(
            db='nuccore', id=nucAcc,
            rettype='fasta', retmode='text',
            strand=strand,
            seq_start=start,
            seq_stop=stop
        ) as seqHandle:
            cds = SeqIO.read(seqHandle, 'fasta')
            if strand == '-':
                cds = cds.reverse_complement()
        cds.name = f'cds_{recAcc}'
        cds.id = f'cds_{recAcc}'
        cds.description = f'{org}_{nucAcc}:{start}-{stop}({strand})'
        resultCdsList.append(cds)
    return resultCdsList


if __name__ == "__main__":
    main()
