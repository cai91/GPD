import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
        
all_uvs_assem={} # uv -> assemblies (non-redundant)
with open('WG_crispr_targets.txt') as inFile: 
    for line in inFile:
        try:
            all_uvs_assem[line.split()[0]].append(line.strip().split()[1])
        except:
            all_uvs_assem[line.split()[0]]=[line.strip().split()[1]]

VCs=[]
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: # No singletons
            VCs.append(toks)
        
X_hq_deep={}
with open('bwa_processed_75_sampleNames.txt') as inFile:
    for line in inFile:
        toks=line.strip().split(',')
        X_hq_deep[toks[0]]=toks[1:]
        
VC_toGenus={}
for idx in range(len(VCs)):
    VC_toGenus[idx]=[]
    
for idx in range(len(VCs)):
    for uv in VCs[idx]:
        assem=all_uvs_assem.get(uv)
        if assem!=None:
            for a in assem:
                VC_toGenus[idx].append(assem_to_genus[a])
    if len(VC_toGenus[idx])!=0:
        VC_toGenus[idx]=list(set(VC_toGenus[idx]))[0]
# I'm mapping uvigs to their VCs
uvs_to_VC={}
for vc_idx in range(len(VCs)):
    for uv in VCs[vc_idx]:
        uvs_to_VC[uv]=vc_idx
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


# S4A

n_conts=[1,2,3,4,5,6]
VCs_conts=[]
VCs_set_glob=[] # Set of VCs found in 1,2...

for my_n in n_conts:
    VCs_glob=[]
    n=0
    with open('VC_continent_span.txt') as inFile:
        for line in inFile:
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

hits=[]
for my_t in ['genus']:
    gVC_genus=assignHost(VCs_set_glob[4]+VCs_set_glob[5],my_t)
    n=0
    z=0
    for vc in gVC_genus.keys():
        if len(gVC_genus[vc])!=0:
            n+=1
        if len(gVC_genus[vc])==1:
            z+=1
    hits.append(z/n)
n=0
z=0
for g in gVC_genus:
    if len(gVC_genus[g])!=0:
        n+=1
        if len(gVC_genus[g])>1: # Is host range able to cross genera?
            z+=1
        
print('Global')
print('Total assigned: '+str(n)+' Host range > 1 genus: '+str(z))

hits=[]
for my_t in ['genus']:
    gVC_genus=assignHost(VCs_set_glob[0],my_t)
    n=0
    z=0
    for vc in gVC_genus.keys():
        if len(gVC_genus[vc])!=0:
            n+=1
        if len(gVC_genus[vc])==1:
            z+=1
    hits.append(z/n)
n=0
z=0
for g in gVC_genus:
    if len(gVC_genus[g])!=0:
        n+=1
        if len(gVC_genus[g])>1: # Is host range able to cross genera?
            z+=1
        
print('Single continent')
print('Total assigned: '+str(n)+' Host range > 1 genus: '+str(z))

sns.set(font_scale=1.2)
sns.set_style("whitegrid")
plt.figure (figsize=(5,5))
plt.ylim(0,0.3)

df_hr=pd.DataFrame()
df_hr['Distribution']=['Single continent VCs','Global VCs']
df_hr['Fraction of VCs with broad host range']=[347/3116.0,36/139.0]
sns.barplot(x='Distribution',y='Fraction of VCs with broad host range',data=df_hr)

plt.show()

# Figure_S4B

# Supplementary figure B

# On average, how many connections each global VC from each order has?
ords=['Bacteroidales','Lachnospirales','Oscillospirales']
order_gs=[]
av_genera=[]

for order_t in ords:
    av_genera.append([])

    vc_osc={}

    with open('Global_dist_hostNet.txt') as inFile:
        for line in inFile:
            vc=line.split(',')[0]
            g=line.strip().split(',')[1]
            if genus_to_order[g]==order_t:
                try:
                    vc_osc[vc].append(g)
                except:
                    vc_osc[vc]=[g]
    for k in vc_osc:
        vc_osc[k]=list(set(vc_osc[k]))

    av=0
    n=0
    avs_bact=[]
    for k in vc_osc:
        av+=len(vc_osc[k])
        avs_bact.append(len(vc_osc[k]))
        av_genera[-1].append(len(vc_osc[k]))
        n+=1

    order_gs.append(av/n)
    
df_orders_gs=pd.DataFrame()
df_orders_gs['Order']=ords
df_orders_gs['Genera per VC']=order_gs
df_orders_gs=df_orders_gs.sort_values(by='Genera per VC',ascending=False)
sns.set(font_scale=1.2)
plt.figure (figsize=(5,5))
sns.set_style("whitegrid")
#plt.ylim(0,4.0)

sns.barplot(x='Order',y='Genera per VC',data=df_orders_gs)

plt.show()