from Bio import SeqIO


input_handle = open("sc04.embl", "rU")
output_handle = open("sc04.gb", "w")


sequences = SeqIO.parse(input_handle, "embl")
count = SeqIO.write(sequences, output_handle, "genbank")


input_handle.close()
output_handle.close()


print("Converted %i records" % count)
