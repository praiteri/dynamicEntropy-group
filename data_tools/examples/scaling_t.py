import re
import numpy as np
import data_tools as qp

data = qp.read_file_to_df("scaling_t.dat")

sizes = [24,96,192,681]
x = []
y = []
l = []
st = []
for s in sizes:
    d = data[data.iloc[:, 0] == s]
    x.append(d.iloc[:, 1])  # Fixed: added [:, 1] to get column
    y.append(d.iloc[:, 2])
    l.append(f"{s}k atoms")
    st.append("p")


for s in sizes:
    d = data[data.iloc[:, 0] == s]
    xx = np.array([np.power(2,n) for n in range(10)])
    x.append(xx)
    y.append(xx * d.iloc[0, 2])
    l.append(None)
    st.append("d")

qp.setup_logger(level='DEBUG')
qp.plot(x, y,
legend=l,
xlabel="# CPUs",
ylabel="ns/day",
title="LAMMPS with Tersoff potential on Setonix",
logx=10,logy=10,
grid=True,
customstyles=st,
colors=[0,1,2,3,0,1,2,3],
#ptype='matplotlib',
ptype='matplotlib',
fout="cpu_perf.png",
linewidth=3,
)
