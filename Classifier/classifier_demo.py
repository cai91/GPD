import numpy as np
import keras
from keras.models import model_from_json

# Load json and create model
json_file = open('model_kgdu.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
model = model_from_json(loaded_model_json)
# Load weights into new model
model.load_weights("model_kgdu.h5")

# Loading features 
X_ph=[]
n=0
with open('KGDU_phages_features.txt') as inFile:
    for line in inFile:
    	# Skipping header
        if n==0:
            n+=1
        else:
            X_ph.append(line.strip().split('\t'))
X_ice=[] 
n=0
with open('KGDU_ICEs_features.txt') as inFile:
    for line in inFile:
    	# Skipping header
        if n==0:
            n+=1
        else:
            X_ice.append(line.strip().split('\t'))

# Prediction phages
out_ph=model.predict_proba(X_ph)
print(out_ph)
# Prediction ices
out_ice=model.predict_proba(X_ice)
print(out_ice)