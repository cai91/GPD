import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency


def assignHost(VCs_query,level):
    # Host range at the species level
    gVC_tax={}
    for vc in VCs_query:
        gVC_tax[vc]=[]

    for i in VCs_query:
        idx=int(i.split('_')[1])
        for uv in VCs[idx]:
            assem=all_uvs_assem.get(uv)
            if assem!=None:
                for a in assem:
                    if level=='species':
                        gVC_tax[i].append(assem_to_spp[a])
                    if level=='genus':
                        gVC_tax[i].append(assem_to_genus[a])
                    if level=='family':
                        gVC_tax[i].append(assem_to_fam[a])
                    if level=='order':
                        gVC_tax[i].append(assem_to_order[a])
                    if level=='class':
                        gVC_tax[i].append(assem_to_class[a])
                    if level=='phylum':
                        gVC_tax[i].append(assem_to_phyla[a])
    for k in gVC_tax.keys():
        gVC_tax[k]=list(set(gVC_tax[k]))
        
    return gVC_tax

# Metadata

scaff_to_gca={}
with open('gca_to_scaf.txt') as inFile:
    for line in inFile:
        scaff_to_gca[line.split()[1].strip()]=line.split()[0]
        
gca_to_scaff={}
for k in list(scaff_to_gca.keys()):
    gca_to_scaff[scaff_to_gca[k]]=k

assem_to_fam={}
assem_to_order={}
assem_to_class={}
assem_to_phyla={}
assem_to_genus={}
assem_to_spp={}

fam_to_phyla={}
order_to_phyla={}
class_to_phyla={}
genus_to_phyla={}
genus_to_fam={}
genus_to_order={}

with open('hgg_bgi_taxonomy.tab') as inFile:
    for line in inFile:
        assem=line.split('\t')[0]
        if len(assem.split('_'))==3:
            assem=assem.split('_')[0]+'_'+assem.split('_')[1]+'#'+assem.split('_')[2]
        elif 'scaffold' in assem:
            assem=scaff_to_gca[assem]
        fam=line.split('\t')[5]
        phyla=line.split('\t')[2]
        order=line.split('\t')[4]
        genus=line.split('\t')[-2]
        classB=line.split('\t')[3]
        spp=line.split('\t')[-1].strip()
        if 'Firmicutes' in phyla:
            phyla='Firmicutes'
            
        assem_to_fam[assem]=fam
        assem_to_order[assem]=order
        assem_to_class[assem]=classB
        assem_to_phyla[assem]=phyla
        assem_to_genus[assem]=genus
        assem_to_spp[assem]=spp
        
        fam_to_phyla[fam]=phyla
        order_to_phyla[order]=phyla
        class_to_phyla[classB]=phyla
        genus_to_phyla[genus]=phyla
        genus_to_fam[genus]=fam
        genus_to_order[genus]=order
        

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

VCs=[]
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: # No singletons
            VCs.append(toks)  

# Assigning uvigs to each sample 
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

# Figure 4A

# crassphage genera
crass_uvs={}
for i in range(1,11):
    crass_uvs['GEN_'+str(i)]=[]
with open('crass_vcs_genera.txt') as inFile:
    for line in inFile:
        if 'GEN_' in line:
            g=line.strip()
        else:
            if 'vig' in line:
                crass_uvs[g].append(line.strip().split()[0]+'_'+line.strip().split()[1])
# Recover uvigs from VCs
for g in crass_uvs:
    for vc in VCs:
        for uv in vc:
            if uv in crass_uvs[g]:
                crass_uvs[g].extend(vc)
                break

clust=[]
for g in crass_uvs:
    clust.append(crass_uvs[g])
    
# crAssphage incidence
conts=['North America','Europe','Oceania','Asia','Africa','South America']
prev_m=[]
for pX_clade in clust:
    prev=[]
    for my_cont in conts:
        NA_n=0
        for k in list(samples_ds_uvs.keys()):
            if bool(set(pX_clade)&set(samples_ds_uvs[k])) and run_toContinent[k]==my_cont:
                NA_n+=1
            n=0
        for k in list(samples_ds_uvs.keys()):
            if run_toContinent[k]==my_cont:
                n+=1
        prev.append(NA_n/n)
    prev_m.append(prev)
    
NA=[]
EU=[]
OC=[]
AS=[]
AF=[]
SA=[]

for i in prev_m:
    NA.append(i[0])
    EU.append(i[1])
    OC.append(i[2])
    AS.append(i[3])
    AF.append(i[4])
    SA.append(i[5])
    
df_stack=pd.DataFrame()
df_stack['North America']=NA
df_stack['Europe']=EU
df_stack['Oceania']=OC
df_stack['Asia']=AS
df_stack['Africa']=AF
df_stack['South America']=SA
df_stack['Genus']=['I','II','III','IV','V','VI','VII','VIII','IX','X']

sns.set(font_scale=1.2)
sns.set_style("whitegrid")
df_stack.set_index('Genus').T.plot(kind='bar', stacked=True,figsize=(10,7))
plt.ylabel('Proportion of samples')
plt.title('Crass-like epidemiology')
plt.xticks(rotation=0)
plt.legend(facecolor='whitesmoke')
plt.show()

# Figure 4B

# Using format VC_idx
uv_toVC={}
VC_toX={}
for v_idx in range(len(VCs)):
    if len(VCs[v_idx])>1:
        for uv in VCs[v_idx]:
            uv_toVC[uv]='VC_'+str(v_idx)
            VC_toX['VC_'+str(v_idx)]=[]
# Samples assigned per VC (usings uvigs as evidence of VC)
with open('bwa_processed_75_sampleNames.txt') as inFile:
    for line in inFile:
        uv=line.split(',')[0]
        toks=line.strip().split(',')[1:] # Samples for each uvig
        my_vc=uv_toVC.get(uv) # Retrieving VC from uvig (VC_idx+1 format)
        if my_vc!=None: # If not a singleton
            VC_toX[my_vc].extend(toks) # Add samples to vc
# Remove redundancy per VC
for vc in VC_toX.keys():
    VC_toX[vc]=list(set(VC_toX[vc]))

conts=['North America','Europe','Asia','Africa','South America','Oceania']

with open('VC_continent_span.txt','w') as outFile:
    outFile.write('VC\tNorth America\tEurope\tAsia\tAfrica\tSouth America\tOceania\n')
    # For each VC check it's global distribution
    for vc in VC_toX.keys():
        outFile.write(vc)
        conts_dir={}
        for c in conts:
            conts_dir[c]=0
        for r in VC_toX[vc]:
            my_c=run_toContinent.get(r)
            if my_c!='NA':
                conts_dir[my_c]+=1
        for c in conts:
            outFile.write('\t'+str(conts_dir[c]))
        outFile.write('\n')

n_conts=[1,2,3,4,5,6]
VCs_conts=[]
VCs_set_glob=[] # Set of VCs found in 1,2...

for my_n in n_conts:
    VCs_glob=[]
    n=0
    with open('VC_continent_span.txt') as inFile:
        for line in inFile:
            # Skip header
            if n==0:
                n+=1
            else:
                z=0
                vc=line.split('\t')[0]
                toks=line.strip().split('\t')[1:]
                for t in toks:
                    if int(t)>0:
                        z+=1
                if z==my_n:    # Change this to control at least or exact
                    VCs_glob.append(vc)
    VCs_conts.append(len(VCs_glob))
    VCs_set_glob.append(VCs_glob)
                
df_contDist=pd.DataFrame()
df_contDist['Continents']=n_conts
df_contDist['Number of VCs']=VCs_conts
df_contDist=df_contDist.sort_values(by='Number of VCs')

# uv -> assemblies (non-redundant)
all_uvs_assem={} 
with open('WG_crispr_targets.txt') as inFile: 
    for line in inFile:
        try:
            all_uvs_assem[line.split()[0]].append(line.strip().split()[1])
        except:
            all_uvs_assem[line.split()[0]]=[line.strip().split()[1]]
VCs_dict={}
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: # No singletons
            VCs_dict[line.strip().split('\t')[0]]=toks 

# Only consider VCs in 5 or more continents
VCs_glob=VCs_set_glob[4]+VCs_set_glob[5]
gVC_genus=assignHost(VCs_glob,'genus')   
gs=[]
for l in gVC_genus.keys():
    for g in gVC_genus[l]:
        gs.append(g)
gs=list(set(gs))

with open('Global_dist_hostNet.txt','w') as outFile:
    for vc in gVC_genus.keys():
        for g in gVC_genus[vc]:
            outFile.write(vc+','+g+'\n')