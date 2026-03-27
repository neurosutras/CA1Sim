"""
Plot Eval #1 features vs targets from recovered_metadata_eval_1_20260322_2311.yaml
Saves both PNG and SVG.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import yaml

# ── Load data (both files have identical features, use 2311)
with open('recovered_metadata_eval_1_20260322_2311.yaml') as f:
    meta = yaml.safe_load(f)

cycles = [1, 2, 3, 4, 5]

sim_areas   = meta['features']['areas']
sim_spikes  = meta['features']['spikes']
sim_troughs = meta['features']['troughs']

tgt_areas   = meta['targets']['areas']
tgt_spikes  = meta['targets']['spikes']
tgt_troughs = meta['targets']['troughs']

total_error = meta['total_error']
log_src     = meta['original_log']
recovered   = meta['recovered_on']

# ── Style
plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        11,
    'axes.spines.top':  False,
    'axes.spines.right':False,
    'axes.linewidth':   1.2,
})

SIM_COLOR  = '#2563EB'   # blue
TGT_COLOR  = '#DC2626'   # red
BG_COLOR   = '#F8FAFC'

x = np.array(cycles)
width = 0.35

fig = plt.figure(figsize=(12, 9), facecolor='white')
fig.suptitle(
    f'Eval #1 — Initial Model vs Target\n'
    f'Total Error: {total_error:,.2f}  |  Source: {log_src}  |  Recovered: {recovered}',
    fontsize=12, fontweight='bold', y=0.98
)

gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.55)

# ── Panel 1: Plateau Area
ax1 = fig.add_subplot(gs[0])
ax1.set_facecolor(BG_COLOR)
b1 = ax1.bar(x - width/2, sim_areas,   width, label='Simulation', color=SIM_COLOR, alpha=0.85, zorder=3)
b2 = ax1.bar(x + width/2, tgt_areas,   width, label='Target',     color=TGT_COLOR, alpha=0.85, zorder=3)
ax1.set_ylabel('Area (mV·s)', fontweight='bold')
ax1.set_title('Plateau Area per Theta Cycle', fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels([f'C{i}' for i in cycles])
ax1.legend(frameon=False)
ax1.yaxis.grid(True, alpha=0.4, zorder=0)
# value labels
for bar in b1:
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
             f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=8.5, color=SIM_COLOR)
for bar in b2:
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
             f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=8.5, color=TGT_COLOR)

# ── Panel 2: Spike Count
ax2 = fig.add_subplot(gs[1])
ax2.set_facecolor(BG_COLOR)
b3 = ax2.bar(x - width/2, sim_spikes,  width, label='Simulation', color=SIM_COLOR, alpha=0.85, zorder=3)
b4 = ax2.bar(x + width/2, tgt_spikes,  width, label='Target',     color=TGT_COLOR, alpha=0.85, zorder=3)
ax2.set_ylabel('Spike Count', fontweight='bold')
ax2.set_title('Spikes per Theta Cycle', fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels([f'C{i}' for i in cycles])
ax2.legend(frameon=False)
ax2.yaxis.grid(True, alpha=0.4, zorder=0)
ax2.set_ylim(0, max(max(sim_spikes), max(tgt_spikes)) + 2)
for bar in b3:
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
             f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=8.5, color=SIM_COLOR)
for bar in b4:
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
             f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=8.5, color=TGT_COLOR)

# ── Panel 3: Inter-burst Trough
ax3 = fig.add_subplot(gs[2])
ax3.set_facecolor(BG_COLOR)
ax3.plot(x, sim_troughs, 'o-',  color=SIM_COLOR, lw=2,   ms=8,  label='Simulation', zorder=3)
ax3.plot(x, tgt_troughs, 's--', color=TGT_COLOR, lw=2,   ms=8,  label='Target',     zorder=3)
ax3.axhline(0, color='gray', lw=0.8, linestyle=':', zorder=2)
ax3.fill_between(x, sim_troughs, tgt_troughs, alpha=0.08, color='purple', zorder=1)
ax3.set_ylabel('Trough Vm (mV)', fontweight='bold')
ax3.set_title('Inter-burst Trough per Theta Cycle', fontweight='bold')
ax3.set_xticks(x)
ax3.set_xticklabels([f'C{i}' for i in cycles])
ax3.legend(frameon=False)
ax3.yaxis.grid(True, alpha=0.4, zorder=0)
for xi, (sv, tv) in enumerate(zip(sim_troughs, tgt_troughs), start=1):
    ax3.annotate(f'{sv:.2f}', (xi, sv), textcoords='offset points', xytext=(-18, 6),
                 fontsize=8, color=SIM_COLOR)
    ax3.annotate(f'{tv:.2f}', (xi, tv), textcoords='offset points', xytext=(4, 6),
                 fontsize=8, color=TGT_COLOR)

# ── Save
for ext in ('svg', 'png'):
    out = f'initial_eval_features_20260322_2311.{ext}'
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f'Saved: {out}')

plt.close()
