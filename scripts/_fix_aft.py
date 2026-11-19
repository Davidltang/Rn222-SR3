import json, sys

# ═══════════════════════════════════════════
# Fix AFT_FV in both correction notebooks
# ═══════════════════════════════════════════

notebooks_to_fix = [
    r'C:/Users/13696/Desktop/Rn222/rn222_s1only-analysis-new.ipynb',
    r'C:/Users/13696/Desktop/Rn222/rn222_s1only-analysis-full-0707.ipynb',
]

for nb_path in notebooks_to_fix:
    with open(nb_path, 'r', encoding='utf-8') as f:
        content = f.read()
    before = content
    content = content.replace('0.13, 0.46', '0.15, 0.455')
    content = content.replace('0.13-0.46', '0.15-0.455')
    content = content.replace('[AFT 0.13-0.455]', '[AFT 0.15-0.455]')  # redundant but safe
    with open(nb_path, 'w', encoding='utf-8') as f:
        f.write(content)
    count = content.count('0.15, 0.455')
    print(f'{nb_path.split("/")[-1]}: {count} replacements')

# Also fix cell_dumps
dump_path = r'C:/Users/13696/Desktop/Rn222/cell_dumps/_cell28.txt'
try:
    with open(dump_path, 'r') as f:
        d = f.read()
    d = d.replace('0.13, 0.46', '0.15, 0.455')
    with open(dump_path, 'w') as f:
        f.write(d)
    print(f'cell_dumps/_cell28.txt: updated')
except FileNotFoundError:
    print('cell_dumps/_cell28.txt: not found, skipping')

print('\nAFT values fixed.')
