#!/usr/bin/env /p/app/python3/bin/python3
# Use this to unpickle the files
import pickle
import json

# Find the best matching model

with open('ranking_debug.json', 'r') as jsonfile:
    data=jsonfile.read()

obj = json.loads(data)

# Unpickle file

filename = 'result_' + obj['order'][0] + '.pkl'
print("unpacking file",filename)
infile = open(filename, 'rb')
new_dict = pickle.load(infile)
infile.close()

#print(new_dict)
print(new_dict['predicted_aligned_error'])
print(len(new_dict['predicted_aligned_error']))

pae=new_dict['predicted_aligned_error']
print(type(pae))

pae_dims=pae.shape
print('pae dim0 = ',pae_dims[0])

with open("pae.csv","w") as csvfile:
    for i in range(pae_dims[0]):
        for j in range(pae_dims[1]):
            csvfile.write("{:.3f} ".format(pae[i,j]))
        csvfile.write("\n")
    
