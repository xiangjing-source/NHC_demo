# ------------------- Simulation control -------------------
units       real          # 单位制: real / lj
dt          1.0           # 时间步长 (fs in real units)
steps       100000        # 总步数
thermo      1           # 每隔多少步输出一次热力学信息
T_target    180.0         # 目标温度

# ------------------- Nose–Hoover chain -------------------
nhc_chain   3             # 链长 M
Tdamp       80.0           # 温控阻尼

# ------------------- Lennard-Jones parameters-------------------
epsilon     0.234         
sigma       3.504         
rc          8.76          
mass        39.948        

# ------------------- File IO -------------------
coords      coords.data   # 初始结构
output      output.log    # 输出日志
plot        md_energy_temp.png   # 输出图片

