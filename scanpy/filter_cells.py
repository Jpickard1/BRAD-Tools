"""
scanpy library: scanpy.pp.filter_cells is used to filter cells during preprocessing

BRAD's Command to Run this file
```
subprocess.run([sys.executable, '<path/to/script/>/filter_cells.py', chatstatus['output-directory'], ], capture_output=True, text=True)
```

Arguments (three arguments):
  1. output directory: chatstatus['output-directory']
  2. output file: <name of output file>
  3. input file: <.h5ad file that can be read to scanpy>
  4. filter type: filter by the maximum or minimum number of genes and cells. Use one of the following: min_counts, max_counts, min_genes, max_genes
  5. filter value: a threshold value for the max/min number of genes or cells to be included in the preprocessed dataset

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `<descriptive name>.pkl`

"""
# Auth: Joshua Pickard
#       jpic@umich.edu
# Date: July 22, 2024

import os
import sys
import scanpy as sc

def main():
    outputPath = sys.argv[1] # chatstatus['output-directory']
    outputFile = sys.argv[2] # output file name
    inputFile  = sys.argv[3] # input file
    filter_type= sys.argv[4] # filter type
    filterVal  = sys.argv[5] # filter value
    filterVal  = int(filterVal)

    adata = sc.read_h5ad(inputFile)
    if filter_type == 'min_counts':
        sc.pp.filter_cells(adata, min_counts=filterVal)
    elif filter_type == 'max_counts':
        sc.pp.filter_cells(adata, max_counts=filterVal)
    elif filter_type == 'min_genes':
        sc.pp.filter_cells(adata, min_genes=filterVal)
    elif filter_type == 'max_genes':
        sc.pp.filter_cells(adata, max_genes=filterVal)
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")

    print('Preprocessing Step ')
    print(' - converted the data with ' + str(filterType) + ' and a threshold of ' + str(filterVal)

    adata.write(os.path.join(outputPath, outputFile))

    print(' - the preprocessed data was saved to ' + os.path.join(outputPath, outputFile))


if __name__ == "__main__":
    main()
