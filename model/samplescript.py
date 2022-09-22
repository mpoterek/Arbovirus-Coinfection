#Marya Poterek
#Script to sample parameter ranges (previously defined) and save in sets of 1K to text files

from SALib.sample import saltelli
import numpy as np

# Define the model inputs
problem = {
    'num_vars': 10,
    'names': ['m', 'a', 'b', 'c', 'r', 'g', 'eiph', 'eipm', 'v2i', 'alpha'],
    'bounds': [[1.0, 4.0],
               [0.3, 1.0],
               [0.05, 0.4],
               [0.25, 0.75],
               [3.0, 10.0],
               [0.05, 0.2],
               [5.0, 8.0],
               [7.0, 11.0],
               [0.0, 75.0],
               [0.0, 1.0]]
}

# Generate samples
param_values = saltelli.sample(problem, 100000)
count = 0
for i in range(0, 2200000, 1000):
    filename = "file" + str(count) + ".txt"
    np.savetxt(filename, param_values[i:i+1000])
    count+=1
