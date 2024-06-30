"""
This script builds linear models of time series data using Dynamic Mode Decomposition (DMD).

Arguments (four arguments):
    1. output directory: chatstatus['output-directory']
    2. output file: <name of output file>
    3. input file: the file created previously in the pipeline with formatted time series data
    4. DMD Rank: this parameter controls the amount reduced order modeling. A value of -1 indicates no model reduction

Usage:
    Command Line:
    ```
    python <path/to/script/>datasetSelector.py <output path> <output file> <input file> <dmd rank>
    ```
                                                     |              |            |           |
                                                 Argument 1     Argument 2   Argument 3  Argument 4
    BRAD Line:
    ```
    subprocess.call([sys.executable, '<path/to/script/>/buildLinearModel.py', chatstatus['output-directory'], <output file>, <input file>, <dmd rank>])
    ```

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `S2-<descriptive name>.pkl`
"""

import pickle
import os
import sys
sys.path.append('/home/jpic/BRAD-Tools-Imports/')
from classes import model # LinearTimeInvariant

def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    inputFile  = sys.argv[3]
    dmdRank    = int(sys.argv[4])
    if dmdRank == -1:
        dmdRank = None
    
    print('******************************')
    print('  Dynamic Mode Decomposition  ')
    print('******************************')
    print(f"Output Path: {outputPath}")
    print(f"Output File: {outputFile}")
    print(f"Input File: {inputFile}")
    print(f"DMD Rank: {dmdRank}")


    # load the input file
    with open(inputFile, 'rb') as file:
        data = pickle.load(file)
    X = data['X']
    states = data['states']
    LTI = model.LinearTimeInvariant(data = X, states = states, dmdRank = dmdRank)
    
    output_file_path = os.path.join(outputPath, outputFile)
    with open(output_file_path, 'wb') as file:
        pickle.dump(LTI, file)
    
    print(f"Model saved to {output_file_path}")

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
