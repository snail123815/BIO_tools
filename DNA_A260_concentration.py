# read seq and a260 form excel and calculate concentration


from DNA_extinction_coefficient import molarConcentration, massConcentration
import pandas as pd
from copy import copy


file = '/Users/durand.dc/Library/Mobile Documents/com~apple~CloudDocs/Resources_iCloud/DNA_real_concentration.xlsx'

df = pd.read_excel(file, index_col='id')

for row in df.iterrows():
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
    molarC = molarConcentration(seq, a260, nType, modifications=modifications)
    massC = massConcentration(seq, a260, nType, modifications=modifications)
    df.loc[name, 'molarC(M)'] = molarC
    df.loc[name, 'massC(g/L)'] = massC
    df.loc[name, 'molarC(μM)'] = molarC * 1000000
    df.loc[name, 'massC(ng/μL)'] = massC * 1000

writer = pd.ExcelWriter(file)
df.to_excel(writer, sheet_name='Sheet1')
for col in writer.sheets['Sheet1'].columns:
    head = col[0].value
    if head.startswith('molarC') or head.startswith('massC'):
        for _cell in col[1:]:
            if _cell.value <= 0.05:
                _cell.number_format = '0.00E+0'
            else:
                _cell.number_format = '0.00'
    if head == 'id':
        for _cell in col[1:]:
            alignment = copy(_cell.alignment)
            alignment.horizontal = 'left'
            _cell.alignment = alignment
writer.save()
