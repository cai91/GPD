import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

# Figure S3A

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
                
X_hq={}
with open('bwa_processed_75_sampleNames.txt') as inFile:
    for line in inFile:
        toks=line.strip().split(',')
        X_hq[toks[0]]=toks[1:]

        
all_samples=run_toStatus.keys()
all_runs_dict={}
# Dictionary with all samples 
for i in all_samples:
    all_runs_dict[i]=0
for i in X_hq.keys():
    for j in X_hq[i]:
        all_runs_dict[j]+=1
        

sns.set(font_scale=1.4)
sns.set_style("whitegrid")
#plt.figure (figsize=(10,5))
plt.figure (figsize=(9,5))
depth=[run_toDepth[i] for i in all_runs_dict.keys()]

plt.scatter(np.array(depth),all_runs_dict.values(),marker='.',color='b',alpha=0.5)

plt.xlabel('Number of reads (x $10^8$)')
plt.ylabel('Number of phages detected')

plt.xticks(np.arange(0,8e8,0.5e8))
plt.yticks(np.arange(0,1001,100))
plt.xlim(0,3e8)
plt.savefig('S3B_1.png',dpi=600,bbox_inches = "tight")
plt.show()