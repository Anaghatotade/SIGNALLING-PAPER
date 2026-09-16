import MDAnalysis as mda
import numpy as np

u = mda.Universe(
    "lammps.data",
    "traj.lammpstrj",
    format="LAMMPSDUMP"
)

print("n atoms =", len(u.atoms))

print("unique resids =", len(np.unique(u.atoms.resids)))
print("first 50 resids =", u.atoms.resids[:50])

print("unique resnums =", len(np.unique(u.atoms.resnums)))
print("first 50 resnums =", u.atoms.resnums[:50])

print("unique segids =", np.unique(u.atoms.segids)[:20])
