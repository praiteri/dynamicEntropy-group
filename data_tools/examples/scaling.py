from data_tools import *
data = read_file_to_df("scaling.dat")

sizes = [30, 82, 240, 660, 810]
x = []
y = []
l = []
for s in sizes:
    d = data[data.iloc[:, 0] == s]
    x.append(d.iloc[:, 1])  # Fixed: added [:, 1] to get column
    y.append(d.iloc[:, 6] / d.iloc[0, 6] / d.iloc[:, 1])
    l.append(f"{s}k atoms")


plot(x, y, 
legend=l, 
xlabel="# GCD",
ylabel="Efficiency",
title="Carbon Ace on Setonix",
ptype='matplotlib',
customstyles=["lp","lp","lp","lp","lp"],
fout="gpu_scaling.png",
)


x = []
y = []
l = []
for s in sizes:
    d = data[data.iloc[:, 0] == s]
    x.append(d.iloc[:, 1])  # Fixed: added [:, 1] to get column
    y.append(d.iloc[:, 2])
    l.append(f"{s}k atoms")


plot(x, y, 
legend=l, 
xlabel="# GCD",
ylabel="Efficiency",
title="Carbon Ace on Setonix",
ptype='matplotlib',
customstyles=["lp","lp","lp","lp","lp"],
fout="gpu_perf.png",
logx=2,logy=10
)
