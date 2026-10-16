import sys,re,collections,csv,bisect
B=''  # run in a folder holding the input files (see README)
chroms={f'HiC_scaffold_{i}' for i in range(1,40)}
# satellite arrays (same rule as the paper: unit 285-305 bp or 2-3x)
sat=collections.defaultdict(list); cur=None
for l in open(B+'B_juncea_trf.txt'):
    if l[0]=='@': cur=l[1:].strip(); continue
    if cur not in chroms: continue
    f=l.split(' ',4); s,e,p=int(f[0]),int(f[1]),int(f[2])
    if any(285*k<=p<=305*k for k in (1,2,3)): sat[cur].append((s,e))
def merge(iv):
    iv.sort(); out=[]
    for s,e in iv:
        if out and s<=out[-1][1]+1: out[-1][1]=max(out[-1][1],e)
        else: out.append([s,e])
    return out
sat={c:merge(v) for c,v in sat.items()}
satbp=sum(e-s+1 for v in sat.values() for s,e in v)
print('satellite bp (chr1-39):',satbp)
# RM .out
rm=collections.defaultdict(list)
for l in open(B+'Belonogaster_juncea.out'):
    f=l.split()
    if len(f)<11 or not f[0].isdigit(): continue
    rm[f[4]].append((int(f[5]),int(f[6]),f[10]))
def cls(c):
    c=c.split('/')[0]
    return {'Simple_repeat':'Simple','Low_complexity':'Low complexity','Satellite':'Satellite','Unknown':'Unclassified'}.get(c,'TE:'+c if c in ('DNA','LTR','LINE','SINE','RC','PLE','Retroposon') else c)
# overlap sat with RM by class (bp, using per-base max over merged intervals per class)
ov=collections.Counter(); anyov=0
for c,iv in sat.items():
    byc=collections.defaultdict(list)
    for s,e,k in rm.get(c,[]): byc[cls(k)].append((s,e))
    allm=merge([(s,e) for s,e,k in rm.get(c,[])])
    def inter(a,b):
        i=j=0;t=0
        while i<len(a) and j<len(b):
            s=max(a[i][0],b[j][0]); e=min(a[i][1],b[j][1])
            if s<=e: t+=e-s+1
            if a[i][1]<b[j][1]: i+=1
            else: j+=1
        return t
    anyov+=inter(iv,allm)
    for k,v in byc.items(): ov[k]+=inter(iv,merge(v))
print('satellite bp masked by RM (any class): %d (%.1f%%)'%(anyov,100*anyov/satbp))
for k,v in ov.most_common(8): print('   ',k,v,'%.1f%%'%(100*v/satbp))
