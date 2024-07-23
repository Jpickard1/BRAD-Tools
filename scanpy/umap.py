"""
Perform Uniform Manifold Approximation and Projection (UMAP) on single-cell RNA-seq data and save the results

BRAD's Command to Run this file:
subprocess.run([sys.executable, '<path/to/script/>/umap.py', chatstatus['output-directory'], <input file>, <output file>, <n_neighbors>, <min_dist>], capture_output=True, text=True)

Arguments:
  1. output directory: chatstatus['output-directory']
  2. input file: <.h5ad file that can be read by scanpy>
  3. output file: <name of output file>
  4. n_neighbors: <number of neighbors for UMAP>
  5. min_dist: <minimum distance parameter for UMAP>

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
    n_neighbors = int(sys.argv[4])  # number of neighbors for UMAP
    min_dist = float(sys.argv[5])  # minimum distance parameter for UMAP

    # Load the input file
    adata = sc.read_h5ad(inputFile)
    
    # Compute UMAP
    sc.tl.umap(adata, n_neighbors=n_neighbors, min_dist=min_dist)
    
    # Save the processed AnnData object to a file
    output_path_full = os.path.join(outputPath, outputFile)
    adata.write(output_path_full)
    print(f'Data written to {output_path_full}')

if __name__ == "__main__":
    main()