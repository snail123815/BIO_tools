
abbDict = {}


with open('journals_edited_DU.txt', 'r') as f:
    for l in f:
        if l == '\n':
            continue
        l = l.strip()
        names = l.split('\t')
        if len(names) == 1: # entry without abbr
            continue
        namesNoQuo = []
        for n in names:
            namesNoQuo.append(n.replace('"',''))
        fullName = namesNoQuo[0].lower()
        abbr = namesNoQuo[1].lower() # only check the first abbr
        if fullName not in abbDict: # entry not found
            abbDict[fullName] = namesNoQuo
            continue
        foundAbbrs = abbDict[fullName]
        if foundAbbrs[0] == abbr: # duplicated
            nfound = len(foundAbbrs)
            nthisl = len(namesNoQuo) -1
            if nfound <= nthisl:
                continue
            else:
                # replace the whole line
                abbDict[fullName] = namesNoQuo
        else:
            if len(foundAbbrs[0]) > len(abbr):
                if foundAbbrs[0].endswith(')'):
                    continue
                else:
                    # replace the whole line
                    abbDict[fullName] = namesNoQuo

with open('journals_mod.txt', 'w') as f:
    for fullName in abbDict:
        l = '\t'.join(abbDict[fullName])
        l +='\n'
        f.write(l)

        
            