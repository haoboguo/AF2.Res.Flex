# Use this to unpickle the files
import pickle
import json

# Find the best matching model

with open('ranking_debug.json', 'r') as jsonfile:
    data=jsonfile.read()

obj = json.loads(data)

# Unpickle file

filename = 'result_' + obj['order'][0] + '.pkl'
infile = open(filename, 'rb')
new_dict = pickle.load(infile)
infile.close()

#print(new_dict)
print(new_dict['plddt'])
print(len(new_dict['plddt']))

with open("plddt.txt", "w") as f:
    for ele in new_dict['plddt']:
        f.write(str(ele) + "\n")
