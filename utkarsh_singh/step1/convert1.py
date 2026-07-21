cell = [
    [12.5150361169283677, 0.0, 0.0],
    [0.0, 12.5150361430482331, 0.0],
    [0.0, 0.0, 12.5150360567636589],
]

with open("CONTCAR.xyz") as fin:
    lines = fin.readlines()

natoms = int(lines[0])

with open("CONTCAR_cart.xyz", "w") as fout:
    fout.write(f"{natoms}\n")
    fout.write("Cartesian coordinates (Angstrom)\n")

    for line in lines[2:]:          # Skip atom count and blank/comment line
        if not line.strip():
            continue

        s, fx, fy, fz = line.split()

        fx = float(fx)
        fy = float(fy)
        fz = float(fz)

        x = fx*cell[0][0] + fy*cell[1][0] + fz*cell[2][0]
        y = fx*cell[0][1] + fy*cell[1][1] + fz*cell[2][1]
        z = fx*cell[0][2] + fy*cell[1][2] + fz*cell[2][2]

        fout.write(f"{s:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
