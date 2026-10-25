import re, io
from collections import Counter
files = ['hawkes_master.tex', 'part1_model_fisher.tex', 'part2_regret.tex', 'part2b_extensions.tex', 'part1b_tangent.tex', 'part3_tracking_review.tex', 'part4_resolutions.tex']
txt = ''.join(io.open(f, encoding='utf8').read() for f in files)
t = re.sub(r'(?<!\\)%.*', '', txt).replace('\\{', '').replace('\\}', '')
depth = 0; bad = None
for i, ch in enumerate(t):
    if ch == '{': depth += 1
    elif ch == '}':
        depth -= 1
        if depth < 0: bad = i; break
print('brace depth at end', depth, '| first negative at', bad)
b = Counter(re.findall(r'\\begin\{(\w+\*?)\}', t)); e = Counter(re.findall(r'\\end\{(\w+\*?)\}', t))
print('environment mismatches:', {k: (b[k], e[k]) for k in set(b) | set(e) if b[k] != e[k]})
labels = Counter(re.findall(r'\\label\{([^}]+)\}', t)); refs = set(re.findall(r'\\ref\{([^}]+)\}', t))
print('duplicate labels:', [k for k, v in labels.items() if v > 1])
print('undefined refs:', sorted(refs - set(labels)))
print('inputs:', re.findall(r'\\input\{([^}]+)\}', t))
