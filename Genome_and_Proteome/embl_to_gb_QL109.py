'''Not only a file converter'''
from Bio import SeqIO
import pathlib  # for checking the existance of a file

# Clear the log file if exist


def clear_log(logfile):
    if pathlib.Path(logfile).is_file():
        log_handle = open(logfile, 'w')
        log_handle.close()
        write_log('Previous log cleaned\n')


def write_log(logstr, end=None):
    log_handle = open(logfile, 'a')
    print(logstr, end=end, file=log_handle)
    print(logstr, end=end)
    log_handle.close()


def embl_to_gb(scnum):
    # Read the file
    embl_file = f"D:/WORKs/Resources/Resource_QL109/Genome_QL109/sc{str(scnum).zfill(2)}.embl"
    output_file = f"D:/WORKs/Resources/Resource_QL109/Genome_QL109/gb_file/sc{str(scnum).zfill(2)}.gb"
    sample = SeqIO.read(embl_file, "embl")
    write_log(f'{embl_file}\n has {len(sample.features)} features.\nConverting to genbank file...\n')

# Parse embl file to gb file with correct annotation:
    write_log('################# Converting features #################\n')
    # write_log('Converting features...\n')
    for feat_num, feat in enumerate(sample.features):
        if 'mol_type' in feat.qualifiers:
            continue

        write_log(f'No.{feat_num}:')

        if feat.qualifiers['gene'][0] != 'unknown':  # if this feature have been assigned to a protein name
            try:
                # Set up 'inference' field to convert convertable information to genbank official format, 'inference' qualifier is more accurate in describing the feature. Here for QL109 we have blastp_match and domain match, so it is a list.
                feat.qualifiers['inference'] = []
                for note in feat.qualifiers['note']:
                    # Every feature has a ORF number in QL109 embl files, write the information on screen, together with protein name and similarity value
                    if note.startswith('ORF '):
                        write_log(note, end=': ')
                        write_log(feat.qualifiers['gene'][0] + '\n\t' +
                                  feat.qualifiers['similarity'][0])
                        ORF_num = note.split(' ')[2]
                        # I will put the information in valid qualifier 'protein_id' and 'locus_tag'
                        feat.qualifiers['note'].remove(note)
                    elif note.startswith('homologous to '):
                        refseq_id = note.split(' ')[2].split('|')[1].split('_')[1]
                        feat.qualifiers['inference'].append(
                            f'similar to AA sequence:UniProtKB:{refseq_id}')
                # Put redundant information in another note
                feat.qualifiers['note'].append(
                    f'blastp_match: {"; ".join(feat.qualifiers["blastp_match"])}, similarity: {"; ".join(feat.qualifiers["similarity"])}')
                del feat.qualifiers["blastp_match"], feat.qualifiers["similarity"]

                if 'domain' in feat.qualifiers:
                    pfams = []
                    for domain in feat.qualifiers['domain']:
                        pfams.append(domain.split(' ')[1])
                        # Put redundant information in another note
                        feat.qualifiers['note'].append(f'domain: {domain}')
                    feat.qualifiers['inference'].append(f'protein motif:PFAM:{",".join(pfams)}')
                    del feat.qualifiers['domain']
            except:
                write_log('Something wrong?')
                raise
        else:  # if the gene has no homologe, calculate orf number here.
            ORF_num = str(feat_num).zfill(4)
            write_log(f'ORF number {ORF_num} on the contig: Unknown protein')

        feat.qualifiers['codon_start'] = 1
        feat.qualifiers['transl_table'] = 11
        locus_tag = f'QL109_{str(scnum).zfill(2)}{ORF_num}'
        protein_id = locus_tag[6:]
        feat.qualifiers['protein_id'] = protein_id
        feat.qualifiers['product'] = feat.qualifiers['gene'][0]
        del feat.qualifiers['gene']
        feat.qualifiers['locus_tag'] = locus_tag
        strand = f'{feat.location.strand:+}'[0]
        write_log(f'\tStrand: {strand}\n')
        # if feat_num == 2:
        # break

    # Rename the sample, if the name is too long, shorten it
    LOCUS = f"QL109_SC{str(scnum).zfill(2)}"
    sample.name = LOCUS

    if pathlib.Path(output_file).is_file():
        caution = 'Caution file overwrote:\n'
    else:
        caution = ''
    write_log(f'Writing genbank file...\n')
    write_flag = SeqIO.write(sample, output_file, "genbank")
    # Show ending information
    if write_flag == 1:
        write_log(f"\n########## Output files ##########\n\n{caution}{output_file}\n{logfile}\n")
    else:
        print(f'Something wrong in writing file to\n{output_file}')


for i in range(1, 14):
    logfile = f"D:/WORKs/Resources/temp/sc{str(i).zfill(2)}_embl_to_gb.log"
    clear_log(logfile)
    embl_to_gb(i)
