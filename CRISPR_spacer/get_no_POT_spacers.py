# POT - potential off target
# If sgRNACas9.pl did not give report, we can parse the POT files to find out those
# spacer without any POT found.
# Result will be saved to npPOT_List.txt

import os

def getSingle(path):
    d = os.path.join(path, 'B.Sort_POT_byID')
    singleHits = []
    for f in os.listdir(d):
        with open(os.path.join(d,f), 'r') as h:
            h.readline() # Skip the first line, it should be the exact hit
            if len(h.readline()) == 0: # if there is no second line, then this one has no POT detected
                singleHits.append(f.split('.')[0])
    singleHits.sort()
    allsp = os.path.join(path, 'report_protospacer_single.txt')
    res = os.path.join(path, 'noPOT_List.txt')
    with open(allsp, 'r') as h:
        with open(res, 'w') as t:
            lines = h.readlines()
            t.write(lines[0])
            for l in lines[1:]:
                spid = l.split('\t')[0].split('.')[1]
                try:
                    singleHits.index(spid)
                except:
                    continue
                t.write(l)