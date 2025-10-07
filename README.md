# NHC_demo

This repository contains a simple molecular dynamics (MD) demo using a Nose–Hoover chain (NHC) thermostat.

Contents
- `md.cpp` : Main C++ MD program (Lennard-Jones particles, periodic box, NHC thermostat).
- `coords.data` : Initial structure in extended XYZ format (atom count and lattice info).
- `in.md` : Text input file with simulation parameters.
- `data.py` : Python script to parse `output.log`, compute statistics, and produce a plot `md_energy_temp.png`.

Build

Requires a C++11 compatible compiler (g++). To compile:

```bash
g++ -std=c++11 -O2 md.cpp -o md
```

Run

1. Edit `in.md` to set parameters (time step, steps, target temperature, etc.).
2. Ensure `coords.data` is present in the same directory.
3. Run the simulation:

```bash
./md in.md
```

This will produce `output.log` and then call `python3 data.py` to generate `md_energy_temp.png`.

Plotting

`data.py` reads `in.md` and `output.log`, prints basic statistics and saves `md_energy_temp.png`.
Matplotlib is required for plotting; install it with:

```bash
pip3 install matplotlib numpy
```

Notes
- The code uses Lennard-Jones interactions with a neighbor list built from a cutoff (`rc`) and a small skin.
- Units: set `units` in `in.md` to `real` or `lj` (the program sets the Boltzmann constant accordingly).
- The input file `in.md` contains human-readable comments; edit values as needed.

License

This repo does not include a license file. Add one if you intend to share or relicense the code.
