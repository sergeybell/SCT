import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = np.loadtxt('WIEN.dat')


#df = pd.DataFrame(data)
#df.columns = ['EB', 'CONS(MU)', 'CONS(NBAR(1))', 'CONS(NBAR(2))', 'CONS(NBAR(3))']
#df.plot(x='EB', y='CONS(MU)')
#df.plot(x='EB', y='CONS(NBAR(1))')
#df.plot(x='EB', y='CONS(NBAR(2))')
#df.plot(x='EB', y='CONS(NBAR(3))')
#plt.show()

def plot_SP(spdat):
    fig, ax = plt.subplots(3,1,sharex=True)
    t = spdat['iteration']
    ax[0].plot(t, spdat['S_X']); ax[0].set_ylabel('S_X')
    ax[1].plot(t, spdat['S_Y']); ax[1].set_ylabel('S_Y')
    ax[2].plot(t, spdat['S_Z']); ax[2].set_ylabel('S_Z')
    ax[2].set_xlabel('turns')
    return fig, ax

TRPRAY = np.loadtxt('TRPRAY.dat')
TRPSPI = np.loadtxt('TRPSPI.dat')

df1 = pd.DataFrame(TRPRAY)
df1.columns = ['iteration', 'ray', 'X', 'A', 'Y', 'B', 'T', 'D']
df2 = pd.DataFrame(TRPSPI)
df2.columns = ['iteration', 'ray', 'S_X', 'S_Y', 'S_Z']
merged_data= df1.merge(df2, on=["iteration","ray"])
print(merged_data)
p1 = merged_data[merged_data['ray'] == 1]
p2 = merged_data[merged_data['ray'] == 2]
p3 = merged_data[merged_data['ray'] == 3]
p4 = merged_data[merged_data['ray'] == 4]
fig, ax = plot_SP(p1)
fig, ax = plot_SP(p2)
fig, ax = plot_SP(p3)
fig, ax = plot_SP(p4)

plt.show()
