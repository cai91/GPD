import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Figure S1C
uv_toCheckV_comp={}

with open('GPD_qc.tsv') as inFile:
    for line in inFile:
        if 'contig_id' not in line:
            toks=line.strip().split('\t')
            uv=toks[0]
            checkV_comp=toks[8]   # checkV_completion
            try:
                uv_toCheckV_comp[uv]=float(checkV_comp)
            except:
                pass
sns.distplot(list(uv_toCheckV_comp.values()))
sns.set_style("whitegrid")
plt.xlabel('Genome completeness (%)')
plt.ylabel('Density')
plt.show()

# Figure S1D
cont_set=[]
with open('GPD_qc.tsv') as inFile:
    for line in inFile:
        if 'contig_id' not in line:
            cont_set.append(float(line.split()[-3]))

sns.set(font_scale=1.3)
sns.set_style("whitegrid")
sns.distplot(cont_set,color='g')

plt.xlabel('Contamination (%)')
plt.ylabel('Density')
plt.show()

# Figure S1E
uvs_clean={}
n=0
with open('all_gpd_cont.tsv') as inFile:
    for line in inFile:
        if 'contig_id' not in line:
            uvs_clean[line.split()[0]]=1

# Size comparison against other databases
size_gpd=[]
size_coll=[]
size_jgi=[]
size_gvd=[]

uvs_sizes={}
with open('sizes_uvigs.fa') as inFile:
    for line in inFile:
        uvs_sizes[line.split()[0]]=int(line.strip().split()[1])
# GPD sizes
for i in list(uvs_clean.keys()):
    size_gpd.append(uvs_sizes[i])
    
with open('FINAL_Gut_Viral_Database_GVD_1.7.2018.fna_sizes.txt') as inFile:
    for line in inFile:
        size_gvd.append(int(line.strip().split()[1]))
with open('JGI_fecal.fa_sizes.txt') as inFile:
    for line in inFile:
        size_jgi.append(int(line.strip().split()[1]))
sns.set(font_scale=1.3)
sns.set_style("whitegrid")
box = plt.boxplot([np.array(size_gpd)/1000.0,np.array(size_jgi)/1000.0,np.array(size_gvd)/1000.0],patch_artist=True,sym='')

colors = ['b', 'lightgrey', 'lightgrey', 'lightgrey']
#colors=['silver']*4
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    
plt.yticks(np.arange(0,101,10))
plt.ylabel('Genome size (kb)')
plt.xlabel('Database')
plt.ylim(0,100)
plt.xticks([1, 2, 3], ['GPD','IMG/VR', 'GVD'])
#plt.title('Genome size distribution')

for item in ['medians']:
    plt.setp(box[item], color='k')
sns.set_style("whitegrid")
plt.show()

# Figure S1F
uvs_phageTax={}
n=0
with open('ivigs_tax-class.tsv') as inFile:
    for line in inFile:
        if n==0:
            n+=1
        else:
            tax=str(line.strip().split('\t'))
            try:
                if float(tax)==1:
                    pass
            except:
                if tax!='':
                    uvs_phageTax[line.split()[0]]=tax
n=0
with open('uvigs_tax-class.tsv') as inFile:
    for line in inFile:
        if n==0:
            n+=1
        else:
            tax=str(line.strip().split('\t'))
            try:
                if float(tax)==1:
                    pass
            except:
                if tax!='':
                    uvs_phageTax[line.split()[0]]=tax
                    
VCs=[]
uv_toVCs={}
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: # No singletons
            VCs.append(toks)
            for t in toks:
                uv_toVCs[t]=line.strip().split('\t')[0]
                
            
fams=['Podoviridae','Siphoviridae','Myoviridae','Herelleviridae','Microviridae','Tectiviridae']
counts_VCs={}
for f in fams:
    counts_VCs[f]=0
for f in fams:
    fam=[]
    for uv in uvs_phageTax:
        if f.lower() in uvs_phageTax.get(uv).lower():
            v=uv_toVCs.get(uv)
            if v!=None:
                fam.append(v)
    fam=set(fam)
    counts_VCs[f]=len(fam)
    
counts_VCs['Unknown']=len(VCs)-sum(counts_VCs.values())
                    
df=pd.DataFrame()
df['Phage family']=counts_VCs.keys()
df['Percentage of VCs']=(np.array(list(counts_VCs.values()))/len(VCs))*100
df=df.sort_values(by='Percentage of VCs',ascending=False)


sns.set(font_scale=1.3)
sns.set_style("whitegrid")
sns.barplot(x='Phage family',y='Percentage of VCs',data=df)
plt.xticks(rotation=90)
plt.yticks(np.arange(0,81,10))
plt.show()