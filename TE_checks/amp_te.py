import csv,collections,re,statistics
B=''  # run in a folder holding the input files (see README)
# gene CDS intervals
cds=collections.defaultdict(list); gchr={}
for l in open(B+'Belonogaster_juncea_final_annot.gff'):
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='CDS': continue
    g=re.search(r'gene_id=([^;]+)',f[8]).group(1)
    cds[g].append((int(f[3]),int(f[4]))); gchr[g]=f[0]
# RM intervals by chrom, TE classes only (interspersed: DNA LTR LINE SINE RC PLE Unknown)
TE=('DNA','LTR','LINE','SINE','RC','PLE','Unknown','Retroposon')
rm=collections.defaultdict(list)
for l in open(B+'Belonogaster_juncea.out'):
    f=l.split()
    if len(f)<11 or not f[0].isdigit(): continue
    c=f[10].split('/')[0]
    if c in TE: rm[f[4]].append((int(f[5]),int(f[6]),f[10]))
for c in rm: rm[c].sort()
starts={c:[s for s,e,k in v] for c,v in rm.items()}
import bisect
def te_frac(g):
    iv=sorted(cds[g]); tot=sum(e-s+1 for s,e in iv); c=gchr[g]; v=rm.get(c,[]); st=starts.get(c,[])
    cov=0; classes=collections.Counter()
    for s,e in iv:
        i=max(0,bisect.bisect_left(st,s)-200)
        mark=set()
        while i<len(v) and v[i][0]<=e:
            a,b,k=v[i]
            if b>=s:
                lo,hi=max(a,s),min(b,e)
                for x in range(lo,hi+1): mark.add(x)
                classes[k.split('/')[0]]+=hi-lo+1
            i+=1
        cov+=len(mark)
    return cov/tot, classes
# families from OrthoFinder v2 N0 (run 15)
rows=list(csv.reader(open(B+'N0_OFv2.tsv'),delimiter='\t'))
h=rows[0]; sp=h[3:]
amp=[];single=[]
for r in rows[1:]:
    cnt={s:len([x for x in (r[3+i] if 3+i<len(r) else '').split(',') if x.strip()]) for i,s in enumerate(sp)}
    bj=[x.strip() for x in (r[3+sp.index('Belonogaster_juncea')]).split(',') if x.strip()]
    o=max(v for k,v in cnt.items() if k!='Belonogaster_juncea')
    if cnt['Belonogaster_juncea']>=20 and cnt['Belonogaster_juncea']>=5*max(o,1): amp.append((r[0],bj))
    if all(v==1 for v in cnt.values()): single+=bj
ampg=[g for _,v in amp for g in v]
print('amplified HOGs:',len(amp),'genes:',len(ampg),'| single-copy genes:',len(single))
def summ(gs,label):
    fr=[];cl=collections.Counter()
    for g in gs:
        if g not in cds: continue
        f,c=te_frac(g); fr.append(f); cl+=c
    n=len(fr)
    print(f'{label}: n={n}; median CDS TE-covered {statistics.median(fr):.2f}; >=50% covered {100*sum(x>=0.5 for x in fr)/n:.0f}%; any overlap {100*sum(x>0 for x in fr)/n:.0f}%; top classes',[(k,round(v/1e3)) for k,v in cl.most_common(5)])
summ(ampg,'amplified'); summ(single,'single-copy')
import random; random.seed(1); summ(random.sample(list(cds),3000),'random 3000 genes')
# per family
for hog,v in sorted(amp,key=lambda x:-len(x[1]))[:12]:
    fr=[te_frac(g)[0] for g in v if g in cds]; cl=collections.Counter()
    for g in v:
        if g in cds: cl+=te_frac(g)[1]
    print(hog,len(v),'median',round(statistics.median(fr),2),'top',cl.most_common(2))
