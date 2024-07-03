"""
This script performs time varying sensor or biomarker selection on dynamical systems models of cellular processes. It selects new sensors at each time point.

Arguments (four arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: the file created previously with the LTI or LTV model
    4. start time: (int, default 0) the first time to select sensors at (inclusive)
    5. end time: (int, default 10) the last time to select sensors at (not inclusive)
    6. verbose: (optional, True or False) if sensor selection should be run in verbose mode so that more output will be displayed to the user.    

Usage:
Command Line:
```
python <path/to/script/>biomarkerSelectionTimeVarying.py <output path> <output file> <input file> <start time> <end time> <verbose>
```

BRAD Line:
```
subprocess.call([sys.executable, '<path/to/script/>/biomarkerSelectionTimeVarying.py', chatstatus['output-directory'], <output file>, <input file>, <start time>, <end time>, <verbose>], capture_output=True, text=True)
```

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `S3-<descriptive name>.csv`
"""

import pickle
import os
import sys
import numpy as np
sys.path.append('/home/jpic/BRAD-Tools-Imports/')
from classes import model # LinearTimeInvariant
import pandas as pd
from scipy.sparse.linalg import eigs

def energyMaximizationTV(MODEL, times, v=False):
    if not isinstance(MODEL, model.LinearTimeInvariant) and not isinstance(MODEL, model.LinearTimeVariant):
        raise ValueError("Unsupported model type. Supported types: LinearTimeInvariant, LinearTimeVariant")
    TVSensors = {}
    for t in times:
        if v:
            print('t: ' + str(t) + '/' + str(len(times)))
        G = MODEL.gram_matrix_TV(T=t)
        # print(G.shape)
        D, V = eigs(G, k=1)
        # D, V = np.linalg.eig(G)
        V = np.abs(V) # this line was missing. it is used in the original Hasnain code.
        obs = pd.DataFrame({'state'  : MODEL.states,
                            'ev1'    : V[:,0],
                            'weight' : np.real(V[:,0])})
        obs['rank'] = obs['weight'].rank(ascending=False)
        obs = obs.sort_values(by='rank', ascending=True)
        obs = obs.reset_index(drop=True)
        TVSensors[t] = {'sensors': obs,
                        'G'      : G,
                        'evals'  : D,
                        'evecs'  : V}
    return TVSensors


def main():
    print('Time varying biomarker selection script has been called!')
    outputPath = sys.argv[1] # str
    outputFile = sys.argv[2] # str
    inputFile  = sys.argv[3] # str
    startTime  = int(sys.argv[4])  # int
    endTime    = int(sys.argv[5])  # int
    if len(sys.argv) > 6:
        verbose = sys.argv[6].lower() in ['true', '1', 't', 'y', 'yes']  # bool
    else:
        verbose = True  # Default value

    print('*******************************')
    print(' Time Varying Sensor Selection ')
    print('*******************************')
    print(f"Output Path: {outputPath}")
    print(f"Output File: {outputFile}")
    print(f"Input File : {inputFile}")
    print(f"startTime  : {startTime}")
    print(f"endTime    : {endTime}")
    print(f"verbose    : {verbose}")

    # Load the model from the input file
    with open(inputFile, 'rb') as file:
        MODEL = pickle.load(file)

    print('..............................')
    print('         Model Loaded         ')

    times = np.arange(startTime, endTime).tolist()
    sensors = energyMaximizationTV(MODEL, times, v=verbose)
    
    print('..............................')
    print('      Sensors Selected        ')
    
    for time in times:
        output_file_path = os.path.join(outputPath, 'time_' + str(time) + "_" + outputFile)
        pd.DataFrame(sensors[time]['sensors']).to_csv(output_file_path)    
        print(f"Model saved to {output_file_path}")

    print('..............................')
    print('      Sensors To File         ')
    

if __name__ == "__main__":
    main()


