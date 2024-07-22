"""
This file orchistrates submitting an array job to the cluster to embed all of the Tabula Sapiens data with scGPT. Run this file when you want to embed all of the Tabula Sapiens dataset.

This file takes no input arguments and it should always be run the exact same way.

Usage:
Command Line:
```
python <path/to/script/>scgpt.py chatstatus['output-directory']
```

BRAD Line:
```
subprocess.run([sys.executable, '<path/to/script/>/scgpt.py', chatstatus['output-directory']], capture_output=True, text=True)
```
"""
# Auth: Joshua Pickard
#       jpic@umich.edu
# Date: July 18, 2024

import subprocess

def run_sbatch():
    # The bash command to run
    command = ['sbatch', '/home/jpic/BRAD-Tools/foundation-models/ts_embedding-scGPT.sh']

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    # Print the output and any errors
    print(f"stdout: {result.stdout}")
    print(f"stderr: {result.stderr}")

if __name__ == "__main__":
    run_sbatch()
