import numpy as np
import sys
from keras.models import model_from_json

file_features=sys.argv[1]

# Load json, create model, load neural network weights
json_file = open('model_kgdu.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
model = model_from_json(loaded_model_json)
model.load_weights("model_kgdu.h5")

# Loading features 
X=[]
with open(file_features) as inFile:
    for line in inFile:
        X.append(line.strip().split('\t'))

# Predicting probability that an input sequence is a phage
out=model.predict_proba(X)
print(out)