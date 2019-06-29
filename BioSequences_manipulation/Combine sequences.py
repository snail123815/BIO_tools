from Bio import SeqIO
import os


def createList(foldername, fulldir=True, suffix=".gb"):
    file_list_tmp = os.listdir(foldername)
    print("Total files", len(file_list_tmp))
    file_list = []
    desired_file_num = 0
    if fulldir:
        for item in file_list_tmp:
            if item.endswith(suffix):
                file_list.append(os.path.join(foldername, item))
                desired_file_num += 1
    else:
        for item in file_list_tmp:
            if item.endswith(suffix):
                file_list.append(item)
                desired_file_num += 1
    print("Files with '" + suffix + "' suffix:", desired_file_num)
    return file_list


file_list = createList(".")


combined = []
for file in file_list:
    seq_obj = SeqIO.read(file, 'gb')
    combined.append(seq_obj)
print(len(combined))


SeqIO.write(combined, 'combined.fa', 'fasta')
