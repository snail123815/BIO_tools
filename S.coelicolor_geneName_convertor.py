import pandas as pd
import numpy as np


strepDB_annotation_excel = '/Users/durand.dc/Documents/works/Resources/Resource_M145/M145_strepDB_annotations.xlsx'

M145_geneTable_tsv = '/Users/durand.dc/Documents/works/Resources/Resource_M145/M145_gene_table.tsv'

uniprot_coelicolor = '/Users/durand.dc/Documents/works/Resources/Resource_M145/S.coelicolor_uniprot.tsv'

input = '/Users/durand.dc/Documents/works/Project_Sporulation/SCO1839_SGR5654/ChIP-Seq_BGI/Peak_calling/Top_genes.txt'

output = '/Users/durand.dc/Documents/works/Project_Sporulation/SCO1839_SGR5654/ChIP-Seq_BGI/Peak_calling/Top_genes_output.txt'

dfUniprot = pd.read_csv(uniprot_coelicolor, delimiter='\t', index_col='Entry')
dfM145GeneTable = pd.read_csv(M145_geneTable_tsv,
                              delimiter='\t', index_col='locus_tag')

listInput = pd.read_csv(input, delimiter='\t')
inputCols = listInput.columns
outputDfs = []
print(inputCols)
for col in inputCols:
    dfOutput = pd.DataFrame(index=listInput.index)
    for idx in listInput.index:
        iId = listInput.loc[idx, col]
        haveQuery = dfUniprot.loc[:, 'Gene names'].str.contains(iId)
        haveQuery.replace(np.nan, False, inplace=True)
        dfResult = dfUniprot[haveQuery]
        foundId = '_'.join(list(dfResult.index))
        dfOutput.loc[idx, f"{col}_uniprot"] = foundId
        if len(dfResult.index) == 0:
            print(iId)
            if iId in dfM145GeneTable.index:
                print(dfM145GeneTable.loc[iId, :])
            else:
                pass
            print('*' * 100)
            print()
        elif len(dfResult.index) > 1:
            print(iId)
            print(dfResult)
            print('*' * 100)
            print()
    dfOutput = pd.concat((listInput, dfOutput), axis=1)
    outputDfs.append(dfOutput)

pd.concat(outputDfs, axis=1).to_csv(output, sep='\t')
