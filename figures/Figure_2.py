import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Figure 2A

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

all_uvs_assem={} 
with open('WG_crispr_targets.txt') as inFile: 
    for line in inFile:
        try:
            all_uvs_assem[line.split()[0]].append(line.strip().split()[1])
        except:
            all_uvs_assem[line.split()[0]]=[line.strip().split()[1]]
            
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
        
# Counting number of assemblies assigned to each VC 
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

# Viral diversity across bacterial clades (genera)

gen_vD={}
gen_set=[]
for i in list(set(genus_to_iso.keys())):
    if i!='NA':
        gen_set.append(i)
gen_set=list(set(gen_set))
        
for f in gen_set:
    gen_vD[f]=0
    
# Counting number of assemblies per genus
for my_gen in gen_set:
    for vc in vc_to_assem.keys():
        n=0
        # For each assembly associated to each VC
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

plt.figure (figsize=(15,5))
sns.set(font_scale=1.2)
sns.set_style("whitegrid")
sns.barplot(x='Genus',y='VCs/isolate',data=genus_dict,hue='Phylum',dodge=False)
plt.xticks(rotation=85)
plt.show()

# Figure 2B

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

def getAssembly_genus(idx,tax):
    for vc in VCs[idx:idx+1]:
        for uv in vc:
            # Preferentially return a prophage assigned assembly
            if 'ivig' in uv:
                assem_iv=iv_to_assem.get(uv)
                if assem_iv!=None:
                    if assem_to_genus.get(assem_iv)==tax:
                            return assem_iv
            # Otherwise use CRISPR data
            assems=all_uvs_assem.get(uv)
            if assems!=None:
                for a in assems:
                    if assem_to_genus.get(a)==tax:
                        return a
        idx+=1

toGCF={}
with open('gcf_mapping.txt') as inFile:
    for line in inFile:
        toGCF[line.split('\t')[0]]=line.split('\t')[1].strip()
toGCA={}
with open('1520_mapping.txt') as inFile:
    for line in inFile:
        toGCA[line.split()[0][1:]]=line.split()[-1].split('.')[0]+'.'+line.split()[-1].split('.')[1]
iv_to_pred={}
with open('all_ivigs.txt') as inFile:
    for line in inFile:
        iv_to_pred[line.split()[0][1:]]=line.strip().split()[1]
iv_to_assem={}
for iv in iv_to_pred.keys():
    iv_to_assem[iv]=contig_assem(iv)

# Cross-genus analysis
VCs_ann=[]
final_conn=[]
cross_genus={}
for idx in range(0,len(VCs)):
    for vc in VCs[idx:idx+1]: # For each VC check host range
        hostR=[]
        for uv in vc:
            assems=all_uvs_assem.get(uv)
            if assems!=None:
                for a in assems:
                    hostR.append(assem_to_genus.get(a))
                    
    hostR=list(set(hostR))
    # If host range is not restricted to single genus
    if len(hostR)>1:
        cross_genus[idx]=hostR
        
for vc in cross_genus:
    assem_1=getAssembly_genus(vc,cross_genus[vc][0])
    assem_2=getAssembly_genus(vc,cross_genus[vc][1])

    if assem_1==None or assem_2==None:
        pass
    else:
        # Make sure pairs belong to same phyla
        if assem_to_phyla[assem_1]==assem_to_phyla[assem_2]:

            # Preparing for iTOL
            if '#' in assem_1:
                assem_1=assem_1.split('#')[0]+'_'+assem_1.split('#')[1]
            if 'GCA' in assem_1:
                assem_1=gca_to_scaff[assem_1]

            if '#' in assem_2:
                assem_2=assem_2.split('#')[0]+'_'+assem_2.split('#')[1]
            if 'GCA' in assem_2:
                assem_2=gca_to_scaff[assem_2]

            final_conn.append(str(assem_1+','+assem_2))
        
# VCs connecting pairs of assemblies in gut isolate tree
final_conn=list(set(final_conn))
for f in final_conn:
    VCs_ann.append(f)
    print(f)