import zipfile, re, shutil
from xml.sax.saxutils import escape
src='hawkes_research_problem_tracker.xlsx'; zin=zipfile.ZipFile(src); x=zin.read('xl/worksheets/sheet1.xml').decode()
add={6:" EXTENSIONS 25 Oct (assumptionM_exact.md, box_cost_and_lower_bound.md, K2_regret_chain.md): (1) Assumption M certified by exact Fourier quadrature for K=1, eps<=0.3, k0 in [0.6,4] (m>=0.032); (2) box cost: polynomial in c+,1/c-, empirically ~s^2 at rho0=1 (prior effect), leading constant box-free; scale of mu* not learnable at eps=c/T; (3) K>=2: chain carries over, Assumption M replaced by computable lambda_L (Thm FK), Thm R^(K) (K+1)/2 ln T + C_K, constants uninformative, K=2 slope ~1.3 vs 1.5 [N]; (4) minimax lower bound (d/2)ln n - O(ln ln n) for ANY predictor, conditional on uniformity (U). LaTeX master document: latex/hawkes_master.tex (not compiled).",
     15:" CLEANUP 25 Oct: stale c_lag / amplitude-independence / K=2 bias-ratio statements struck in the note; noise convention (sigma^2=a=r/(2 eff)); path-length bound V^2/(eta T) (factor 1/2 not derived); jump term units. Label clashes resolved in LABELS.md."}
for row,t in add.items():
    m=re.search(rf'<c r="M{row}"( s="\d+")[^>]*?(?:/>|><is><t[^>]*>(.*?)</t></is></c>)',x,re.S); assert m,row
    x=x[:m.start()]+f'<c r="M{row}"{m.group(1)} t="inlineStr"><is><t>{m.group(2) or ""}{escape(t)}</t></is></c>'+x[m.end():]
tmp='tracker_new.xlsx'
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zo:
    for i in zin.infolist(): zo.writestr(i, x.encode() if i.filename=='xl/worksheets/sheet1.xml' else zin.read(i.filename))
zin.close(); shutil.move(tmp,src); print(zipfile.ZipFile(src).testzip())
