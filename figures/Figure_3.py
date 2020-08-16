import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing.data import StandardScaler
from sklearn.decomposition import PCA

# Figure 3A

# Assign VCs in format VC_idx 
vc_to_uvs={}
vc_idx=0
VCs=[]
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: 
            VCs.append(toks)
            vc_to_uvs[vc_idx]=toks
            vc_idx+=1

uv_to_VC={}
for k in vc_to_uvs:
    uvs=vc_to_uvs[k]
    for uv in uvs:
        uv_to_VC[uv]=k
        
# Mapping samples to VCs
X_to_VCs={}
with open('DeepSamples_x_uvFormat.txt') as inFile:
    for line in inFile:
        toks=line.strip().split(',')
        l=[]
        for t in toks[1:]: # t is each uvig 
            v=uv_to_VC.get(t) # Convert each mapped uvig to an existing VC (VC_idx + 1 format)
            if v!=None:
                l.append(v)
        l=list(set(l)) # Remove redundancy
        X_to_VCs[toks[0]]=l

# Metadata
scaff_to_gca={}
with open('gca_to_scaf.txt') as inFile:
    for line in inFile:
        scaff_to_gca[line.split()[1].strip()]=line.split()[0]
        
gca_to_scaff={}
for k in list(scaff_to_gca.keys()):
    gca_to_scaff[scaff_to_gca[k]]=k
        
all_uvs_assem={} 
with open('WG_crispr_targets.txt') as inFile: 
    for line in inFile:
        try:
            all_uvs_assem[line.split()[0]].append(line.strip().split()[1])
        except:
            all_uvs_assem[line.split()[0]]=[line.strip().split()[1]]


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

# Figure 3A

# Only samples with deep sequencing taken into account for analysis
all_samples=run_toStatus.keys()
samples_ds=[]
for i in all_samples:
    if run_toDepth[i]>=0.5e8:
        v=X_to_VCs.get(i)
        if v!=None:
            samples_ds.append(i)

# Analysing only deep sequenced samples
h_run=[]
for x in samples_ds:
    if run_toStatus[x]=='Healthy':
        if run_toLife[x]=='Urban' or run_toLife[x]=='Rural':
            if run_toContinent[x]=='Asia' or run_toContinent[x]=='Africa' or run_toContinent[x]=='South America' or run_toContinent[x]=='Europe' or run_toContinent[x]=='North America' or run_toContinent[x]=='Oceania': 
                h_run.append(x)

X_all=h_run # List of samples to analyze

with open('diff_X.txt','w') as outFile:
    for k1 in X_all:
        for k2 in X_all:
            vc1=set(X_to_VCs[k1]) # This maps to VC_idx + 1 
            vc2=set(X_to_VCs[k2])
            try:
                outFile.write(k1+'\t'+k2+'\t'+str(len(vc1&vc2)/float((len(vc1)+len(vc2))))+'\n')
            except:
                outFile.write(k1+'\t'+k2+'\t'+str(0)+'\n')

X_d=[]
n=len(X_all)
with open('diff_X.txt') as inFile:
    for line in inFile:
        if n==len(X_all):
            n=0
            X_d.append([])
            X_d[-1].append(float(line.split('\t')[-1]))
            n+=1
        else:
            X_d[-1].append(float(line.split('\t')[-1]))
            n+=1

plt.figure (figsize=(7,5))

scaler = StandardScaler()
X_scaled = scaler.fit(X_d).transform(X_d)

pca = PCA(n_components=2)
X_r = pca.fit(X_scaled).transform(X_scaled)

X_rx=[i[0] for i in X_r]
X_ry=[i[1] for i in X_r]


country_sp=[]
for x in h_run:
    if run_toCountry[x]=='Fiji' or run_toCountry[x]=='United Republic of Tanzania' or run_toCountry[x]=='Madagascar' or run_toCountry[x]=='Peru':
        if run_toCountry[x]=='United Republic of Tanzania':
            country_sp.append('Tanzania')
        else:
            country_sp.append(run_toCountry[x])
    else:
        country_sp.append('Other')
        
lifestyle_sp=[]
for x in h_run:
    if run_toCountry[x]=='Fiji' or run_toCountry[x]=='United Republic of Tanzania' or run_toCountry[x]=='Madagascar' or run_toCountry[x]=='Peru':
        lifestyle_sp.append('Rural')
    else:
        lifestyle_sp.append(run_toLife[x])

df_pca=pd.DataFrame()
df_pca['PC1 (18.24%)']=X_rx
df_pca['PC2 (5.35%)']=X_ry
df_pca['Run']=h_run
df_pca['Status']=len(h_run)*['Healthy']
df_pca['Study']=[run_toPub[x] for x in h_run]
df_pca['Continent']=[run_toContinent[x] for x in h_run]
df_pca['Country']=country_sp
df_pca['Age']=[run_toAge[x] for x in h_run]
df_pca['Sample lifestyle']=lifestyle_sp

plt.figure (figsize=(9,7))
sns.set(font_scale=1.3)
sns.scatterplot(x='PC1 (18.24%)',y='PC2 (5.35%)',data=df_pca,s=60,hue='Sample lifestyle',alpha=0.65,palette=['b','g','grey'],hue_order=['Urban','Rural','NA']) #style='Continent'
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

# Figure 3B

VC_toGenus={}
for idx in range(len(VCs)):
    VC_toGenus[idx]=[]
for idx in range(len(VCs)):
    for uv in VCs[idx]:
        assem=all_uvs_assem.get(uv)
        if assem!=None:
            for a in assem:
                VC_toGenus[idx].append(assem_to_genus[a])
    VC_toGenus[idx]=list(set(VC_toGenus[idx]))
    
# Assigning uvigs to each deep sample
X_hq_deep={}
with open('bwa_processed_75_sampleNames.txt') as inFile:
    for line in inFile:
        toks=line.strip().split(',')
        X_hq_deep[toks[0]]=toks[1:]
all_samples=run_toStatus.keys()
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

# Bacteroides targeting

conts_prev=[]
sns.set(font_scale=1.0)
conts=['North America','Europe','Asia','Africa','South America']

for c in conts:
    genera=[]
    S=[]
    for s in samples_ds:
        if run_toContinent[s]==c:
            for uv in samples_ds_uvs[s]:
                my_vc=uv_to_VC.get(uv)
                if my_vc!=None:
                    g=VC_toGenus.get(my_vc)
                    genera.extend(g)

    g_dict={}
    gen_set=list(set(genera))
    for i in gen_set:
        g_dict[i]=0
    for i in genera:
        g_dict[i]+=1
        
    all_s=sum(g_dict.values())
    for i in g_dict:
        g_dict[i]=g_dict[i]/all_s
        
    #ts=['Bacteroides','Bacteroides_A','Bacteroides_B','Bacteroides_F']
    ts=['Bacteroides','Bacteroides_A','Bacteroides_B']
        
    s_prop=0
    for t in ts:
        v=g_dict.get(t)
        if v!=None:
            s_prop+=v
    
    conts_prev.append(s_prop)


country_prev=[]
sns.set(font_scale=1.0)
country=['Australia','Fiji']

for c in country:
    genera=[]
    S=[]
    for s in samples_ds:
        if run_toCountry[s]==c:
            for uv in samples_ds_uvs[s]:
                my_vc=uv_to_VC.get(uv)
                if my_vc!=None:
                    g=VC_toGenus.get(my_vc)
                    genera.extend(g)
                    
                    
    g_dict={}
    gen_set=list(set(genera))
    for i in gen_set:
        g_dict[i]=0
    for i in genera:
        g_dict[i]+=1
        
    all_s=sum(g_dict.values())
    for i in g_dict:
        g_dict[i]=g_dict[i]/all_s
        
    #ts=['Bacteroides','Bacteroides_A','Bacteroides_B','Bacteroides_F']
    ts=['Bacteroides','Bacteroides_A','Bacteroides_B']
        
    s_prop=0
    for t in ts:
        v=g_dict.get(t)
        if v!=None:
            s_prop+=v
    
    country_prev.append(s_prop)

sns.set(font_scale=1.5)
plt.figure (figsize=(9,5))

df_prev=pd.DataFrame()
df_prev['Population']=conts[0:3]+country[0:1]+conts[3:]+country[1:]
df_prev['Proportion of viral sequences']=conts_prev[0:3]+country_prev[0:1]+conts_prev[3:]+country_prev[1:]

sns.set_style("whitegrid")
sns.barplot(x='Population',y='Proportion of viral sequences',data=df_prev,color='b')
plt.yticks(np.arange(0,0.61,0.1))
plt.title('Bacteroides')
plt.xticks(rotation=80)
plt.show()

print(df_prev['Proportion of viral sequences'])

# Prevotellaceae

conts_prev=[]
sns.set(font_scale=1.0)
conts=['North America','Europe','Asia','Africa','South America']

for c in conts:
    #print(c)
    genera=[]
    S=[]
    for s in samples_ds:
        if run_toContinent[s]==c:
            for uv in samples_ds_uvs[s]:
                my_vc=uv_to_VC.get(uv)
                if my_vc!=None:
                    g=VC_toGenus.get(my_vc)
                    genera.extend(g)
                    #if g!=None:
                        #if g!=[]:
                            #genera.append(g)
    g_dict={}
    gen_set=list(set(genera))
    for i in gen_set:
        g_dict[i]=0
    for i in genera:
        g_dict[i]+=1
        
    all_s=sum(g_dict.values())
    for i in g_dict:
        g_dict[i]=g_dict[i]/all_s
        
    ts=['Prevotella','Paraprevotella']
        
    s_prop=0
    for t in ts:
        v=g_dict.get(t)
        if v!=None:
            s_prop+=v
    
    conts_prev.append(s_prop)


# **Gut virome (bacterial host)**
country_prev=[]
sns.set(font_scale=1.0)
country=['Australia','Fiji']
#prev_p=[]
for c in country:
    #print(c)
    genera=[]
    S=[]
    for s in samples_ds:
        if run_toCountry[s]==c:
            for uv in samples_ds_uvs[s]:
                my_vc=uv_to_VC.get(uv)
                if my_vc!=None:
                    g=VC_toGenus.get(my_vc)
                    genera.extend(g)
    g_dict={}
    gen_set=list(set(genera))
    for i in gen_set:
        g_dict[i]=0
    for i in genera:
        g_dict[i]+=1
        
    all_s=sum(g_dict.values())
    for i in g_dict:
        g_dict[i]=g_dict[i]/all_s
        
    ts=['Prevotella','Paraprevotella']
        
    s_prop=0
    for t in ts:
        v=g_dict.get(t)
        if v!=None:
            s_prop+=v
    
    country_prev.append(s_prop)

sns.set(font_scale=1.5)
plt.figure (figsize=(9,5))

df_prev=pd.DataFrame()
df_prev['Population']=conts[0:3]+country[0:1]+conts[3:]+country[1:]
df_prev['Proportion of viral sequences']=conts_prev[0:3]+country_prev[0:1]+conts_prev[3:]+country_prev[1:]

sns.set_style("whitegrid")
sns.barplot(x='Population',y='Proportion of viral sequences',data=df_prev,color='g')
plt.yticks(np.arange(0,0.61,0.1))
plt.title('Prevotellaceae')
plt.xticks(rotation=80)
plt.show()

print(df_prev['Proportion of viral sequences'])