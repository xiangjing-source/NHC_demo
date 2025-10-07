# ------------------- Simulation control -------------------
units       real          # units: real / lj
dt          1.0           # time step (fs in real units)
steps       100000        # total number of steps
thermo      1              # output thermodynamic info every this many steps
T_target    180.0         # target temperature

# ------------------- Nose–Hoover chain -------------------
nhc_chain   3             # chain length M
Tdamp       80.0          # thermostat damping time

# ------------------- Lennard-Jones parameters-------------------
epsilon     0.234         
sigma       3.504         
rc          8.76          
mass        39.948        

# ------------------- File IO -------------------
coords      coords.data   # initial structure file
output      output.log    # output log file
plot        md_energy_temp.png   # output image file

