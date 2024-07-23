"""
Perform differential expression analysis on single-cell RNA-seq data and save the results

BRAD's Command to Run this file:
subprocess.run([sys.executable, '<path/to/script/>/differential_expression.py', chatstatus['output-directory'], <input file>, <output file>, <groupby>], capture_output=True, text=True)

Arguments:
  1. output directory: chatstatus['output-directory']
  2. input file: <.h5ad file that can be read by scanpy>
  3. output file: <name of output file>
  4. groupby: <column name in adata.obs to group by>

**OUTPUT FILE NAME INSTRUCTIONS**
1. Output path should be chatstatus['output-directory']
2. Output file name should be `<descriptive name>.h5ad`
"""

# Auth: Joshua Pickard
#       jpic@umich.edu
# Date: July 22, 2024

import os
import sys
import scanpy as sc

def main():
    outputPath = sys.argv[1]  # chatstatus['output-directory']
    inputFile = sys.argv[2]   # input file
    outputFile = sys.argv[3]  # output file name
    groupby = sys.argv[4]    # column name in adata.obs to group by

    # Load the input file
    adata = sc.read_h5ad(inputFile)
    
    # Perform differential expression analysis
    sc.tl.rank_genes_groups(adata, groupby=groupby)
    
    # Save the processed AnnData object to a file
    output_path_full = os.path.join(outputPath, outputFile)
    adata.write(output_path_full)
    print(f'Data written to {output_path_full}')

if __name__ == "__main__":
    main()
