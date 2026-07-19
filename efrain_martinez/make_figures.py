from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


materials = ["Sn-Br", "Sn-I", "Ti-Br", "Ti-I", "Zr-Br", "Zr-I"]
m0_s1 = np.array([5.216844, 3.368605, 4.942565, 3.715114, 6.237475, 4.809327])
m1_s1 = np.array([0.329880, 0.155191, 0.203623, 0.030523, 0.091242, 0.374627])
m1_bright = np.array([0.531794, 0.209239, 0.203623, 0.420172, 0.091242, 0.541928])
m1_bright_f = np.array([0.001451, 0.077379, 0.053920, 0.004628, 0.011504, 0.193408])

out = Path("Figures")
out.mkdir(exist_ok=True)

x = np.arange(len(materials))
w = 0.36

fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.bar(x - w / 2, m0_s1, w, label="M0")
ax.bar(x + w / 2, m1_s1, w, label="M1")
ax.set_ylabel("State-1 excitation energy (eV)")
ax.set_xticks(x, materials)
ax.set_ylim(0, 6.8)
ax.legend(frameon=False)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig(out / "figure_1_M0_M1_state1_energy.png")
plt.close(fig)

fig, ax1 = plt.subplots(figsize=(7.2, 4.2))
ax1.bar(x - w / 2, m1_s1, w, label="M1")
ax1.bar(x + w / 2, m1_bright, w, label="Brightest")
ax1.set_ylabel("Excitation energy (eV)")
ax1.set_xticks(x, materials)
ax1.set_ylim(0, 0.62)
ax1.grid(axis="y", alpha=0.25)

ax2 = ax1.twinx()
ax2.plot(x, m1_bright_f, "o-", label="Brightest f", color="black")
ax2.set_ylabel("Oscillator strength")
ax2.set_ylim(0, 0.22)
handles = ax1.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
labels = ax1.get_legend_handles_labels()[1] + ax2.get_legend_handles_labels()[1]
ax1.legend(handles, labels, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(out / "figure_2_M1_state1_vs_brightest.png")
plt.close(fig)

print("Wrote Figures/figure_1_M0_M1_state1_energy.png")
print("Wrote Figures/figure_2_M1_state1_vs_brightest.png")
