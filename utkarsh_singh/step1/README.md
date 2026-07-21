# Step 1: Run molecular dynamics using CP2K

## What is this

 - 'CONTCAR` - geometry and unit cell information, usually obtained from VASP cell-relaxation run

 - `submit.slm` - SLURM submit file


## How to run

1. Obtain 'CONTCAR` from VASP cell-relaxation

2. Obtain `.xyz` file and unit cell information:


This will generate the Cartesian coordinates in `TiO2_unit_cell.xyz` and will print out the unit cell informaiton. Use it to prepare MD input file.

3. Edit `md.inp` file as needed

4. Edit `submit.slm` file as needed

5. Run the calculations on the HPC:

```bash
sbatch submit.slm
```

## Results

- `*-pos-1.xyz` file contains the trajectory - WE WILL NEED it in the following step

- `*-vel-1.xyz` velocities

- `*-frc-1.xyz` forces 

- `*.log` file - general output of CP2



Note: This part was computed on my JNCASR's HPC facility using VASP 6.2.0. But, one can do it on UB CCR following the above instructions.K

