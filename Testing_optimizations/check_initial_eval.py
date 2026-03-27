import matplotlib.pyplot as plt
import numpy as np
import optimize_soma_plateaus as osp

print("\n--- Running Single Initial Evaluation ---")
# Use the starting parameters
error = osp.objective(osp.x_start)

print(f"\nInitial Error: {error:.4f}")

# Access the best_trace from the module to see the update
if osp.best_trace is not None:
    time_axis = np.arange(0, osp.duration, osp.dt)
    plt.figure(figsize=(12, 6))
    
    # Slice trace to match time axis if needed
    trace = osp.best_trace[:len(time_axis)]
    
    plt.plot(time_axis, trace, color='black', linewidth=1.5, label='Initial Soma Trace')
    plt.axhline(y=-70, color='gray', linestyle='--', alpha=0.5)
    plt.axhline(y=0, color='red', linestyle='--', alpha=0.3)
    
    # Add target indicators for AUC windows
    for start in [500, 650, 800, 950, 1100]:
        plt.axvspan(start, start + 150, color='blue', alpha=0.05)

    plt.title(f'Initial Model State (Error: {error:.2f})')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.legend()
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    
    save_path = 'initial_eval_trace_fixed.png'
    plt.savefig(save_path, dpi=300)
    plt.savefig('initial_eval_trace_fixed.svg')
    print(f"Plot saved as: {save_path} and initial_eval_trace_fixed.svg")
    # plt.show() # Can't show in terminal
else:
    print("Error: No trace recorded.")
