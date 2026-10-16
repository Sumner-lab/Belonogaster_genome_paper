"""Re-test the between-species overlap of caste-biased genes at the orthogroup level
(first done with orthogroup_analysis.2.R) with two universes:
  A) all orthogroups in the orthogroup file
  B) orthogroups containing >=1 tested gene (non-NA padj in DESeq2) in BOTH species of the pair (stricter)
DEGs: padj < 0.05; LFC>0 = queen-biased (reproductive), LFC<0 = worker-biased."""
import csv,math,itertools,sys
OGFILE=sys.argv[1]; OUT=sys.argv[2]
from scipy.stats import hypergeom
base=''  # folder holding the DESeq2 result tables
res={'BJ':base+'res_Bjun 1.tsv','LF':base+'res_LF.tsv','PD':base+'res_Pdom.tsv'}
col={'BJ':'Belonogaster_juncea','LF':'Liostenogaster_flavolineata','PD':'Polistes_dominula'}
tested={};deg={}
for sp,f in res.items():
    t=set();q=set();w=set()
    for r in csv.DictReader(open(f),delimiter='\t'):
        if r['padj'] in ('NA',''): continue
        t.add(r['gene'])
        if float(r['padj'])<0.05:
            (q if float(r['log2FoldChange'])>0 else w).add(r['gene'])
    tested[sp]=t; deg[(sp,'Queen')]=q; deg[(sp,'Worker')]=w
    print(sp,'tested',len(t),'queen-biased',len(q),'worker-biased',len(w))
og={}
for r in csv.DictReader(open(OGFILE),delimiter='\t'):
    og[r['Orthogroup']]={sp:{g.strip().split('.t')[0] for g in r[c].split(',') if g.strip()} for sp,c in col.items()}
print('orthogroups',len(og))
def ogs(sp,genes): return {o for o,d in og.items() if d[sp]&genes}
OG={k:ogs(k[0],v) for k,v in deg.items()}
TOG={sp:ogs(sp,tested[sp]) for sp in col}
rows=[]
pairs=[(a,b) for a,b in itertools.combinations(['BJ','LF','PD'],2)]
tests=[((a,c),(b,c)) for a,b in pairs for c in ('Queen','Worker')]+[((a,'Worker'),(b,'Queen')) for a,b in pairs]+[((a,'Queen'),(b,'Worker')) for a,b in pairs]
for (k1,k2) in tests:
    for uname,U in (('A_all',set(og)),('B_tested_both',TOG[k1[0]]&TOG[k2[0]])):
        s1=OG[k1]&U; s2=OG[k2]&U; ov=len(s1&s2); N=len(U)
        exp=len(s1)*len(s2)/N; p=hypergeom.sf(ov-1,N,len(s1),len(s2))
        rows.append([f'{k1[0]} {k1[1]}',f'{k2[0]} {k2[1]}',uname,N,len(s1),len(s2),ov,round(exp,1),round(ov/exp,2) if exp else '',p])
# BH within each universe for the 6 same-caste tests
for uname in ('A_all','B_tested_both'):
    sc=[r for r in rows if r[2]==uname and r[0].split()[1]==r[1].split()[1]]
    ps=sorted(range(len(sc)),key=lambda i:sc[i][-1]); m=len(sc); prev=1
    adj=[0]*m
    for rank,i in reversed(list(enumerate(ps,1))):
        prev=min(prev,sc[i][-1]*m/rank); adj[i]=prev
    for r,a in zip(sc,adj): r.append(a)
w=csv.writer(open(OUT,'w'),delimiter='\t')
w.writerow(['set1','set2','universe','N','size1','size2','overlap','expected','fold','P','BH_same_caste'])
for r in rows: w.writerow(r)
for r in rows:
    print('\t'.join(str(x if not isinstance(x,float) else ('%.3g'%x)) for x in r))
