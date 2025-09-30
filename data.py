import numpy as np
import matplotlib.pyplot as plt

# ------------------- 从 in.md 读取参数 -------------------
def read_params(filename="in.md"):
    params = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # 去掉行内注释
            line = line.split("#")[0].strip()
            parts = line.split()
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1]
                try:
                    params[key] = float(val)
                except ValueError:
                    pass  # 如果不是数字就忽略
    return params

params = read_params("in.md")
T_target = params.get("T_target", 300.0)  # 默认 300K

# ------------------- 读取模拟结果 -------------------
data = np.loadtxt("output.log", comments="#", skiprows=1)

steps = data[:, 0]
epot = data[:, 1]
ekin = data[:, 2]
etot = data[:, 3]
temp = data[:, 4]

# ------------------- 统计 -------------------
avg_temp = np.mean(temp)
rms_temp = np.sqrt(np.mean((temp - T_target)**2))
avg_ekin = np.mean(ekin)
rms_ekin = np.sqrt(np.mean((ekin - np.mean(ekin))**2))

print(f"Average temperature: {avg_temp:.5f}")
print(f"Temperature RMS deviation from target: {rms_temp:.5f}")
print(f"Average kinetic energy: {avg_ekin:.5f}")
print(f"Kinetic energy RMS deviation: {rms_ekin:.5f}")

# ------------------- 绘图 -------------------
fig, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

# 动能曲线
axes[0].plot(steps, ekin, label="E_kin")
axes[0].set_ylabel("Energy")
axes[0].set_title("Energy vs Steps")
axes[0].legend(loc='upper right')  # 图例右上角
axes[0].grid(True)
# 标注平均值和 RMS（左上角）
axes[0].text(0.02, 0.95, f"Avg E_kin = {avg_ekin:.5f}\nRMS = {rms_ekin:.5f}",
             transform=axes[0].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

# 温度曲线
axes[1].plot(steps, temp, label="Temperature", color="tab:red")
axes[1].axhline(T_target, color="tab:blue", linestyle="--", label=f"Target Temp = {T_target}")
axes[1].set_xlabel("Step")
axes[1].set_ylabel("Temperature")
axes[1].set_title("Temperature vs Steps")
axes[1].legend(loc='upper right')  # 图例右上角
axes[1].grid(True)
# 标注平均温度和 RMS（左上角）
axes[1].text(0.02, 0.95, f"Avg Temp = {avg_temp:.5f} \nRMS = {rms_temp:.5f} ",
             transform=axes[1].transAxes,
             verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6))

plt.tight_layout()
plt.savefig("md_energy_temp.png", dpi=300)
