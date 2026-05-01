
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nested.optimize_utils import OptimizationHistory
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python3 plot_optimization.py <hdf5_file>")
    sys.exit(1)

hdf5_path = sys.argv[1]
if not os.path.exists(hdf5_path):
    print(f"File not found: {hdf5_path}")
    sys.exit(1)

print(f"Loading optimization history from {hdf5_path}...")
history = OptimizationHistory(file_path=hdf5_path)

print("Generating plots...")
# We plot objectives to see convergence and parameters to see distributions
# For a small test run, we'll just plot the main objective scores
history.plot(subset=['objectives'])

fignums = plt.get_fignums()
print(f"Found {len(fignums)} figures. Saving...")

for i in fignums:
    plt.figure(i)
    # The first few figures are usually: 1. Fitness Rank, 2. Multi-objective score, 3. Total Error
    # Then one for each objective/parameter
    title = plt.gca().get_title()
    filename = f"optimization_results_{i}.png"
    if title:
        safe_title = "".join([c if c.isalnum() else "_" for c in title])
        filename = f"optimization_{safe_title}.png"
    
    plt.savefig(filename, bbox_inches='tight')
    print(f"Saved {filename}")

plt.close('all')
print("Done.")
