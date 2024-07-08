"""
This script organizes datasets into a format to work with time series gene expression using individual gene expresession as the state space. From a .h5ad file, a .pkl file is constructed here that will contain the state type being gene expression, the state names being gene names, and the time series data organized as states (genes) by timepoint (hours) by replicates (experiment replicates).

Arguments (three arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: <file created in previous step>

Usage:

Command Line:
```
python <path/to/script/>geneCoordinateTimeSeries.py <output path> <output file> <input file>
```
                                                         |              |           |
                                                     Argument 1     Argument 2   Argument 3
BRAD Line:
```
subprocess.run([sys.executable, '<path/to/script/>/geneCoordinateTimeSeries.py', chatstatus['output-directory'], <output file>, <input file>], capture_output=True, text=True)
```

*Always replace <path/to/script> with the correct path given above.*

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `S1-<data set name>.pkl`
"""
import os  # All BRAD tools
import sys # All BRAD tools

import anndata as an
import numpy as np
import pickle


def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    dataset    = sys.argv[3] # input file

    print('******************************')
    print('  Convert to Gene Coordinates ')
    print('******************************')
    
    print('Dataset=', end='')
    print(dataset)

    ad = an.read(dataset)
    
    # Build the object in which to save the trajectories. This will be state variables by time by trajectory
    n = ad.shape[1]
    T = ad.obs['order'].nunique()
    R = ad.obs['replicate'].nunique()
    X = np.zeros((n, T, R))

    # List possible timepoints & replicates
    timepoints = ad.obs['order'].unique()
    replicates = ad.obs['replicate'].unique()

    # Loop over dataframe
    for i in range(ad.obs.shape[0]):
        order = ad.obs['order'][i]
        replicate = ad.obs['replicate'][i]
    
        time = np.where(timepoints == order)[0][0]
        rep = np.where(replicates == replicate)[0][0]
    
        X[:, time, rep] = ad.X[i, :]
    states = list(ad.var.index)
    metadata = ad.var
    
    output_file = os.path.join(outputPath, outputFile)
    with open(output_file, 'wb') as f:
        pickle.dump({'X': X, 'states': states, 'metadata': metadata}, f)
    
    print(f"Data saved to {output_file}")


if __name__ == "__main__":
    main()