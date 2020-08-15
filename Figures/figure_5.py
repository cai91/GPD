import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Figure 5A
uvs_gpd={}
with open('GPD_summary.tsv') as inFile:
    for line in inFile:
        if 'contig_id' not in line:
            uvs_gpd[line.split()[0]]=1

sns.set(font_scale=1.2)
VCs=[]
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        l=[]
        toks=line.strip().split('\t')
        for t in toks:
            v=uvs_gpd.get(t)
            if v!=None:
                l.append(t)
        if len(l)>1: # No singletons
            VCs.append(l)
            
sort_VCs=[len(i) for i in VCs]
sort_VCs.sort(reverse=True)
        
f, ax = plt.subplots(figsize=(8, 4))
y=sort_VCs[0:100]
l=list(range(1,101))
x_s=['1','','','','','','','','','10','','','','','','','','','','20','','','','','','','','','','30','','','','','','','','','','40','','','','','','','','','','50','','','','','','','','','','60','','','','','','','','','','70','','','','','','','','','','80','','','','','','','','','','90','','','','','','','','','','100']
ax.plot(l,y,c='grey',alpha=0.7)
ax.scatter(l,y,marker='.')

ax.set_xticks(l)
ax.set_xticklabels(x_s)

plt.yticks(np.arange(0,2000,100))
plt.xlabel('Viral cluster')
plt.ylabel('Genomes/VC')
plt.grid(axis='x')
plt.scatter([1],sort_VCs[0],c='r',marker='*',s=80.0)
plt.scatter([2],sort_VCs[1],c='r',marker='*',s=80.0)


# Figure 5B
X=[]
X_ids=[]
with open('pX_matrix.txt') as inFile:
    for line in inFile:
        uv=line.split(',')[0]
        toks=[float(i) for i in line.split(',')[1:]]

        X.append(toks)
        X_ids.append(uv)

            
# pX structure
sns.clustermap(X,vmin=0.4,vmax=1.2,cmap=plt.cm.Blues)
plt.show()
plt.show()

# Figure 5C

# Fetching metadata
n=1
run_toCountry={}
run_toContinent={}
run_toStatus={}
run_toDepth={}
run_toDisease={}
run_toPub={}
run_toAge={}
run_toStudy={}
run_toLife={}
with open('Gut-metagenomes_29052019.csv') as inFile:
    for line in inFile:
        if n==1:
            n+=1
        else:
            if line.split(',')[2]=='Yes':
                my_run=line.split(',')[0]
                run_toCountry[my_run]=line.split(',')[13]
                run_toContinent[my_run]=line.split(',')[14]
                run_toStatus[my_run]=line.split(',')[5]
                run_toDepth[my_run]=float(line.split(',')[1])
                run_toDisease[my_run]=line.split(',')[6]
                run_toPub[my_run]=line.strip().split(',')[-1]
                run_toLife[my_run]=line.strip().split(',')[12]
                
                age=line.split(',')[9]
                if age!='NA':
                    run_toAge[my_run]=float(age)
                else:
                    run_toAge[my_run]='NA'
                run_toStudy[my_run]=line.split(',')[4]
                
all_samples=run_toStatus.keys()

X_hq_deep={}
with open('bwa_processed_75_sampleNames.txt') as inFile:
    for line in inFile:
        toks=line.strip().split(',')
        X_hq_deep[toks[0]]=toks[1:]

samples_ds=[]
for i in all_samples:
    if run_toDepth[i]>=0.5e8:
        samples_ds.append(i)

samples_ds_uvs={}
for s in samples_ds:
    samples_ds_uvs[s]=[]
    
for k in list(X_hq_deep.keys()):
    for s in X_hq_deep[k]:
        samples_ds_uvs[s].append(k)
        
all_pX_uvs=[]
with open('pX.txt') as inFile:
    for line in inFile:
        all_pX_uvs.append(line.strip())

conts=['North America','Europe','Oceania','Asia','Africa','South America']

cont_pX={}
for c in conts:
    cont_pX[c]=0
    
for c in conts:
    n=0
    for s in samples_ds_uvs:
        if run_toContinent[s]==c:
            n+=1
            if set(samples_ds_uvs[s])&set(all_pX_uvs):
                cont_pX[c]+=1
    cont_pX[c]=cont_pX[c]/float(n)
    
df_pX=pd.DataFrame()
df_pX['Continent']=cont_pX.keys()
df_pX['Prevalence (%)']=np.array(list(cont_pX.values()))*100
df_pX=df_pX.sort_values(by='Prevalence (%)',ascending=False)

sns.set(font_scale=1.4)
sns.barplot(x='Continent',y='Prevalence (%)',data=df_pX,color='b')#palette='Reds_r')
plt.yticks(np.arange(0,41,5))
#plt.title('shortPodo family prevalence around the world')
plt.xticks(rotation=85)
plt.show()
