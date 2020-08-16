import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def contig_assem(iv):
    '''This function links a contig with its assembly'''
    b=1
    h=iv_to_pred[iv]
    if h.split('_')[1]=='' or h.split('_')[1][0:2]=='ER':
        my_id=h.split('_')[2]+'_'+h.split('_')[3]+'_'+h.split('_')[4]
        assem=my_id.split('_')[0]+'_'+my_id.split('_')[1]+'#'+my_id.split('_')[2]
        b=0
    if h.split('_')[1]=='NZ' or h.split('_')[1]=='NC':
        my_id=h.split('_')[1]+'_'+h.split('_')[2]
        assem=toGCF[my_id]
        b=0
    if b:
        try:
            sc=scaff_to_gca[h.split('_Scaf')[0].split('_')[-1]+'_scaffold']
            assem=sc
        except:
            print(iv)
    return assem

def assignHost(VCs_query,level):
    # Host range at the species level for a list of VCs
    gVC_tax={}
    for vc in VCs_query:
        gVC_tax[vc]=[]

    for i in VCs_query: # For each query VC
        idx=int(i.split('_')[1])
        for uv in VCs[idx]: # Examine each uvig from that VC
            assem=all_uvs_assem.get(uv) # Assign assembly to each uvig
            if assem!=None:
                for a in assem: # Assign taxa at the requested rank
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
                        
    # For each query VC set it to a non-redundant list
    for k in gVC_tax.keys():
        gVC_tax[k]=list(set(gVC_tax[k]))
        
    return gVC_tax

# Figure S2A
VCs=[]
with open('GPD_VCs.txt') as inFile:
    for line in inFile:
        toks=line.strip().split('\t')
        if len(toks)>1: # No singletons
            VCs.append(toks)  

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

# Number of VCs / number of isolates (phyla)
phyla_to_iso={}
phyla_set=list(set(assem_to_phyla.values()))
for k in phyla_set:
    phyla_to_iso[k]=0
    
for k in phyla_set:
    with open('hgg_bgi_taxonomy.tab') as inFile:
        for line in inFile:
            k_i=line.split('\t')[2]  # Phyla
            if 'Firmicutes' in k_i:  # This is to simplify Firmicutes_X to just Firmicutes
                k_i='Firmicutes'
            if k_i==k:
                phyla_to_iso[k]+=1

all_uvs_assem={} # uv -> assemblies (non-redundant)
with open('WG_crispr_targets.txt') as inFile: 
    for line in inFile:
        try:
            all_uvs_assem[line.split()[0]].append(line.strip().split()[1])
        except:
            all_uvs_assem[line.split()[0]]=[line.strip().split()[1]]

# Infected isolates per phyla

phyla_set=list(set(assem_to_phyla.values()))
phyla_count={}
for p in phyla_set:
    phyla_count[p]=[]
    
for uv in list(all_uvs_assem.keys()):
    for a in all_uvs_assem[uv]:
        p=assem_to_phyla[a]
        if p!=None:
            phyla_count[p].append(a)
            
for p in phyla_set:
    phyla_count[p]=(len(set(phyla_count[p]))/phyla_to_iso[p])*100.0
    
# Removal of keys due to low counts

del phyla_count['Desulfobacterota']
del phyla_count['Synergistota']
del phyla_count['Campylobacterota']
del phyla_count['Fusobacteriota']

phyla_dict=pd.DataFrame()
phyla_dict['Phyla']=phyla_count.keys()
phyla_dict['Isolates linked to phage (%)']=phyla_count.values()
phyla_dict=phyla_dict.sort_values(by='Isolates linked to phage (%)',ascending=False)

sns.set(font_scale=1.4)
sns.set_style("whitegrid")
sns.barplot(x='Phyla',y='Isolates linked to phage (%)',data=phyla_dict,color='m')
plt.yticks(np.arange(0,101,10))
plt.xticks(rotation=85)
plt.show()

# Figure S2B
# Viral diversity by phylum

# Counting number of predictions assigned to each phyla (viral diversity)
vc_to_assem={}
for vc_idx in range(len(VCs)):
    vc_to_assem[vc_idx]=[]
    for uv in VCs[vc_idx]:
        assem=all_uvs_assem.get(uv)
        if assem!=None:
            for assem_i in assem: 
                vc_to_assem[vc_idx].append(assem_i)        
for k in vc_to_assem.keys():
    vc_to_assem[k]=list(set(vc_to_assem[k]))

phyla=['Firmicutes','Proteobacteria','Bacteroidota','Actinobacteriota']
phyla_vd={}
for p in phyla:
    phyla_vd[p]=[]
    
for vc in vc_to_assem:
    for a in vc_to_assem[vc]:
        v=phyla_vd.get(assem_to_phyla[a])
        if v!=None:
            phyla_vd[assem_to_phyla[a]].append(vc)
for p in phyla_vd:
    phyla_vd[p]=list(set(phyla_vd[p]))
    phyla_vd[p]=len(phyla_vd[p])/phyla_to_iso[p]
    
phyla_dict=pd.DataFrame()
phyla_dict['Phyla']=phyla_vd.keys()
phyla_dict['VCs/isolate']=phyla_vd.values()
phyla_dict=phyla_dict.sort_values(by='VCs/isolate',ascending=False)

sns.set(font_scale=1.4)
sns.set_style("whitegrid")
sns.barplot(x='Phyla',y='VCs/isolate',data=phyla_dict,color='b')
plt.yticks(np.arange(0,3.6,0.5))
plt.xticks(rotation=85)
plt.show()


# Figure S2C
  
# Number of VCs / number of isolates (genus)
genus_to_iso={}
genus_set=list(set(assem_to_genus.values()))
for k in genus_set:
    genus_to_iso[k]=0

for k in genus_set:
    with open('hgg_bgi_taxonomy.tab') as inFile:
        for line in inFile:
            if line.split('\t')[-2]==k:
                genus_to_iso[k]+=1  

sns.set(font_scale=1.4)
my_c=1
VCs_glob=['VC_'+str(i) for i in range(len(VCs))]
hits=[]
for my_t in ['species','genus','family','order','class','phylum']:
    gVC_tax=assignHost(VCs_glob,my_t) # Get taxa at the requested rank for all VCs
    n=0
    z=0
    for vc in gVC_tax.keys(): # For each VC
        if len(gVC_tax[vc])!=0: # This is the host assignment for the VC
            n+=1
        if len(gVC_tax[vc])==1: # How many VCs are contained within 1 taxonomic rank
            z+=1
    hits.append(z/n)

hits=[0]+hits
hits_d=[]
for h_idx in range(1,len(hits)):
    hits_d.append(hits[h_idx]-hits[h_idx-1])

X=hits_d

df_hr=pd.DataFrame()
df_hr['Phylogenetic barrier']=['Species','Genus','Family','Order','Class','Phylum']
df_hr['Percentage of VCs']=np.array([X[0],X[1],X[2],X[3],X[4],X[5]])*100.0
sns.set_style("whitegrid")
sns.barplot(x='Phylogenetic barrier',y='Percentage of VCs',data=df_hr,palette='Blues_d')
plt.yticks(np.arange(0,71,5))
plt.show()

df_hr

# Figure S2D
# Genus level analysis
n_i=0
tot=0
lc=[]
super_spp=[]
for idx in range(0,len(VCs)):
    for vc in VCs[idx:idx+1]: # For each VC check host range
        hostR=[]
        for uv in vc:
            assems=all_uvs_assem.get(uv)
            if assems!=None:
                for a in assems:
                    hostR.append(assem_to_genus.get(a))
    hostR=list(set(hostR))
    # If host range is not restricted to genus
    if len(hostR)>0:
        super_spp.append(hostR)


gen_vD={}
gen_set=[]
for i in list(set(genus_to_iso.keys())):
    if i!='NA':
        gen_set.append(i)
gen_set=list(set(gen_set))
        
for f in gen_set:
    gen_vD[f]=0
    
for my_gen in gen_set:
    for vc in vc_to_assem.keys():
        n=0
        for assem in vc_to_assem[vc]:
            gen_i=assem_to_genus.get(assem)
            if gen_i==my_gen:
                n+=1
        gen_vD[my_gen]+=n
    gen_vD[my_gen]=gen_vD[my_gen]/genus_to_iso[my_gen]
    
# Only picking families with at least 10 members
gen_vD_n10={}
for g in gen_vD.keys():
    if genus_to_iso[g]>=10:
        gen_vD_n10[g]=gen_vD[g]
# Adding phylum
phylum=[]
for g in gen_vD_n10.keys():
    phylum.append(genus_to_phyla[g])
        
genus_dict=pd.DataFrame()
genus_dict['Genus']=gen_vD_n10.keys()
genus_dict['VCs/isolate']=gen_vD_n10.values()
genus_dict['Phylum']=phylum
genus_dict=genus_dict.sort_values(by='VCs/isolate',ascending=False)

all_gen_set=gen_vD_n10.keys()

genTo_bhr={}
for g in all_gen_set:
    genTo_bhr[g]=0
    n=0
    for i in super_spp:
        if (g in i) and len(i)>1:
            n+=1
    genTo_bhr[g]+=n
    
sns.set_style("whitegrid")
plt.scatter([x for x in gen_vD_n10.values()],[y for y in genTo_bhr.values()],marker='.',color='k')
plt.xlabel('Viral diversity')
plt.ylabel('Broad host range hits')
plt.show()