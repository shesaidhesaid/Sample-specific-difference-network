import scipy.stats as stat
import numpy as np
#读入背景网络
ref_1 = {}
ref_m = []
s = []
f = open("E:\\Bio\\workrectify\\27342\\gene_27342_raw_normal.txt")
flag=0
n=0
for p in f:
    flag += 1
    t = p.split()
    if flag==1:
        n=len(t)
        s=t[1:]
        continue
    ref_1[t[0]]=[float(t[i]) for i in range(1,n)]
    ref_m.append([float(t[i]) for i in range(1,n)])
f.close()
genes = []
for i in ref_1.keys():
    genes.append(i)

ref_m=np.matrix(ref_m)
ref_m=np.corrcoef(ref_m)
#读入样本
sample = {}
f=open("E:\\Bio\\workrectify\\27342\\gene_27342_raw.txt")
flag=0
n=0
for p in f:
    flag += 1
    t = p.split()
    if flag==1:
        n=len(t)
        name_d=t[1:]
        continue
    sample[t[0]]=[float(t[i]) for i in range(1,n)]
f.close()

def ssn_score(deta,pcc,nn):
    if pcc == 1:
        pcc = 0.99999999
    if pcc == -1:
        pcc = -0.99999999
    z = deta/((1-pcc*pcc)/(nn-1))
    return z

for k in range(20,50):
    sample_m = []
    fw = open("E:\\Bio\\workrectify\\27342\\ref_27_27_0.01\\"+str(k)+"_27_value.txt","w")
    for s in range(len(genes)):
        sample_m.append(ref_1[genes[s]]+[sample[genes[s]][k]])
    sample_m=np.matrix(sample_m)
    sample_m=np.corrcoef(sample_m)
    for i in range(len(genes)-1):
        for j in range(i+1,len(genes)):
            r = sample_m[i,j]-ref_m[i,j]
            z = ssn_score(r,ref_m[i,j],len(ref_1[genes[i]]))
            pvalue = 1-stat.norm.cdf(abs(z))
            if float(pvalue)<0.01:
                fw.write(genes[i]+"\t"+genes[j]+"\n")
    fw.close()





    
