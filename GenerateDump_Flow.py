import numpy as np
import os
import os.path
import sys

N = int(sys.argv[1])
mu_c = int(sys.argv[2])
mu_e = int(sys.argv[3])
v0_c = int(sys.argv[4])
v0_e = int(sys.argv[5])
pid = int(sys.argv[6])

mu_str = "Mu_" + str(mu_c) + "E" + str(mu_e)
v0_str = "v0_" + str(v0_c) + "E" + str(v0_e)
pid_str = str(pid).zfill(5)

filedir = "/Users/wangd/Documents/Yale/FrictionVibrationModes/RiverBed/Programs/Flow/Files/N_" + str(N).zfill(4) + "/" + mu_str + "/" + pid_str
##filedir = "/Users/wangd/Documents/Yale/DPM/CPP/Adipocyte/3D/"
flowdir = filedir + "/" + v0_str + "/"
sedname = filedir + "/Pos_Sediment_" + pid_str + ".txt"
dumpname = filedir + "/PosFlow_All_" + mu_str + "_" + pid_str + "_" + v0_str + ".dump"

Dl = 1.4
Ds = 1.0
D_old = np.zeros(N)
Nl = int(N / 2)
for n in range(Nl):
    D_old[n] = Dl
    D_old[n + Nl] = Ds

pos_sed = np.loadtxt(sedname)
L = pos_sed[0]
x_old = pos_sed[1:2*N+1:2]
y_old = pos_sed[2:2*N+1:2]

N_up = 0
N_bottom = 0
x_new = np.zeros(N)
y_new = np.zeros(N)
D_new = np.zeros(N)
for n in range(N):
    if y_old[n] < 0.5 * D_old[n]:
        N_bottom += 1
        j = N - N_bottom
        x_new[j] = x_old[n]
        y_new[j] = y_old[n]
        D_new[j] = D_old[n]
    else:
        x_new[N_up] = x_old[n]
        y_new[N_up] = y_old[n]
        D_new[N_up] = D_old[n]
        N_up += 1
R_new = D_new / 2.0

## Nstep = len([f for f in os.listdir(flowdir) if f.startswith('Pos_Flow_') and os.path.isfile(os.path.join(flowdir, f))])
posname = flowdir + "/Pos_Flow_All.txt"
pos_all = np.loadtxt(posname)
pos_count = 2 * N_up
Nstep = int(pos_all.shape[0] / pos_count)

f_dump = open(dumpname, 'w')

for it in range(Nstep):
    pos = pos_all[it*pos_count:(it+1)*pos_count:]
    x = pos[0::2]
    y = pos[1::2]

    for n in range(N_up):
        x[n] -= np.floor(x[n] / L) * L

    f_dump.write("ITEM: TIMESTEP\n")
    f_dump.write(f"{it:d}\n")
    f_dump.write("ITEM: NUMBER OF ATOMS\n")
    f_dump.write(f"{N:d} \n")
    f_dump.write("ITEM: BOX BOUNDS pp pp fixed\n")
    f_dump.write(f"{0:.16e} {L:.16e}\n")
    f_dump.write(f"{0:.16e} {L:.16e}\n")
    f_dump.write(f"{-1:.16e} {1:.16e}\n")
    f_dump.write("ITEM: ATOMS id type mol radius xu yu zu\n")

    for n in range(N_up):
        f_dump.write(f"{n+1:d} 1 {n+1:d} {R_new[n]:.12f} {x[n]:.12f} {y[n]:.12f} {0:.12f}\n")

    for n in range(N_up, N):
        f_dump.write(f"{n+1:d} 1 {2*N_up:d} {R_new[n]:.12f} {x_new[n]:.12f} {y_new[n]:.12f} {0:.12f}\n")

f_dump.close()
