import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Figure 1B

# Custom colors
custom_p=[(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
 (0.8666666666666667, 0.5176470588235295, 0.3215686274509804),
 (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
 (0.5058823529411764, 0.4470588235294118, 0.7019607843137254),
 (0.5764705882352941, 0.47058823529411764, 0.3764705882352941),
 (0.8549019607843137, 0.5450980392156862, 0.7647058823529411),
 (0.5490196078431373, 0.5490196078431373, 0.5490196078431373),
 (0.8, 0.7254901960784313, 0.4549019607843137),
 (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]

# Dictionary to collect checkV hits according to category
qc_cats={}
for i in ['Complete','High-quality','Medium-quality','Low-quality','Not-determined']:
    qc_cats[i]=0
    
# Collecting data from checkV
with open('checkV_GPD_summary.tsv') as inFile:
    for line in inFile:
        if 'checkv_quality' not in line:
            uv=line.split()[0]
            qc=line.split()[6]
            qc_cats[qc]+=1

# Putting data into dataframe
df_qc=pd.DataFrame()
df_qc['High quality']=[qc_cats['Complete'],qc_cats['High-quality'],0,0,0]
df_qc['Genome fragment']=[0,0,qc_cats['Medium-quality'],qc_cats['Low-quality'],qc_cats['Not-determined']]
df_qc['CheckV quality']=['Complete','High quality','Medium quality','Low quality','Not-determined']

# Plotting
sns.set(font_scale=1.2)
sns.set_palette(custom_p)
df_qc.set_index('CheckV quality').T.plot(kind='bar', stacked=True,figsize=(5,7))
plt.ylabel('Number of predictions')
plt.xticks(rotation=0)
plt.yticks(np.arange(0,105001,10000))
plt.show()