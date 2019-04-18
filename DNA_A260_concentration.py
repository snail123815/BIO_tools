# read seq and a260 form excel and calculate concentration


from DNA_extinction_coefficient import molarConcentration, massConcentration
import pandas as pd
from copy import copy


file = '/Users/durand.dc/Library/Mobile Documents/com~apple~CloudDocs/Resources_iCloud/DNA_real_concentration.xlsx'

df = pd.read_excel(file, index_col='id')

molarCM = pd.Series(index=df.index, name='molarC(M)')
massCgl = pd.Series(index=df.index, name='massC(g/L)')
molarCuM = pd.Series(index=df.index, name='molarC(μM)')
massCngul = pd.Series(index=df.index, name='massC(ng/μL)')
for i, row in enumerate(df.iterrows()):
    name = row[0]
    seq, nType, a260 = (row[1]['seq'], row[1]['type'], row[1]['A260'])
    modifications = row[1]['modifications']
    if type(modifications) == str:
        modifications = modifications.split(', ')
    else:
        modifications = [None, ]
    if pd.isnull(a260):
        continue
    print(f"\n\n{'*'*60}\n{name}")
    # print(a260)
    molarC = molarConcentration(seq, a260, nType, modifications=modifications)
    massC = massConcentration(seq, a260, nType, modifications=modifications)
    molarCM.iloc[i] = molarC
    massCgl.iloc[i] = massC
    molarCuM.iloc[i] = molarC * 1000000
    massCngul.iloc[i] = massC * 1000

resultCols = ['molarC(M)', 'massC(g/L)', 'molarC(μM)', 'massC(ng/μL)']
oriDataCols = [col for col in df.columns if col not in resultCols]
newDf = pd.concat((df.loc[:, oriDataCols],
                   molarCM, massCgl, molarCuM, massCngul),
                  axis=1)


writer = pd.ExcelWriter(file)
newDf.to_excel(writer, sheet_name='Sheet1')
# excel formatting
for col in writer.sheets['Sheet1'].columns:
    head = col[0].value
    if head.startswith('molarC') or head.startswith('massC'):
        for _cell in col[1:]:
            if type(_cell.value) != float:
                continue
            if _cell.value <= 0.05:
                _cell.number_format = '0.00E+0'
            else:
                _cell.number_format = '0.00'
    if head == 'id':
        for _cell in col[1:]:
            alignment = copy(_cell.alignment)
            alignment.horizontal = 'left'
            _cell.alignment = alignment
for row in writer.sheets['Sheet1'].rows:
    for _cell in row:
        alignment = copy(_cell.alignment)
        alignment.horizontal = 'left'
        _cell.alignment = alignment
    break
writer.save()
