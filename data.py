import numpy as np
import matplotlib.pyplot as plt

# ------------------- Read parameters from in.md -------------------
def read_params(filename="in.md"):
    params = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # remove inline comments
            line = line.split("#")[0].strip()
            parts = line.split()
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1]
                try:
                    params[key] = float(val)
                except ValueError:
                    pass  # ignore non-numeric values
    return params

params = read_params("in.md")
T_target = params.get("T_target", 300.0)  # default 300K

# ------------------- Read simulation results -------------------
data = np.loadtxt("output.log", comments="#", skiprows=1)

steps = data[:, 0]
epot = data[:, 1]
ekin = data[:, 2]
etot = data[:, 3]
temp = data[:, 4]

# ------------------- Statistics -------------------
avg_temp = np.mean(temp)
rms_temp = np.sqrt(np.mean((temp - T_target)**2))
avg_ekin = np.mean(ekin)
rms_ekin = np.sqrt(np.mean((ekin - np.mean(ekin))**2))

print(f"Average temperature: {avg_temp:.5f}")
print(f"Temperature RMS deviation from target: {rms_temp:.5f}")
print(f"Average kinetic energy: {avg_ekin:.5f}")
print(f"Kinetic energy RMS deviation: {rms_ekin:.5f}")

# ------------------- Plotting -------------------
fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

# Kinetic energy curve
axes[0].plot(steps, ekin, label="E_kin")
axes[0].set_ylabel("Energy")
axes[0].set_title("Energy vs Steps")
axes[0].legend(loc='upper right')  # legend at upper right
axes[0].grid(True)
# annotate average and RMS (upper-left)
axes[0].text(0.02, 0.95, f"Avg E_kin = {avg_ekin:.5f}\nRMS = {rms_ekin:.5f}",
             transform=axes[0].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

# Temperature curve
axes[1].plot(steps, temp, label="Temperature", color="tab:red")
axes[1].axhline(T_target, color="tab:blue", linestyle="--", label=f"Target Temp = {T_target}")
axes[1].set_xlabel("Step")
axes[1].set_ylabel("Temperature")
axes[1].set_title("Temperature vs Steps")
axes[1].legend(loc='upper right')  # legend at upper right
axes[1].grid(True)
# annotate average temperature and RMS (upper-left)
axes[1].text(0.02, 0.95, f"Avg Temp = {avg_temp:.5f} \nRMS = {rms_temp:.5f} ",
             transform=axes[1].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

plt.tight_layout()
plt.savefig("md_energy_temp.png", dpi=300)
