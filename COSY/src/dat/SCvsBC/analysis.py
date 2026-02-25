import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_3(data, X, Y1, Y2, Y3):
    fig, ax = plt.subplots(3,1,sharex=True)
    t = data[X]
    ax[0].plot(t, data[Y1]); ax[0].set_ylabel(Y1)
    ax[1].plot(t, data[Y2]); ax[1].set_ylabel(Y2)
    ax[2].plot(t, data[Y3]); ax[2].set_ylabel(Y3)
    ax[2].set_xlabel(X)
    return fig, ax

def plot_all(fringe_field):

    parameters_X = pd.read_table('Reference/X/Params_'+fringe_field+'.dat', sep="\s+")
    parameters_Y = pd.read_table('Reference/Y/Params_'+fringe_field+'.dat', sep="\s+")
    parameters_D = pd.read_table('Reference/D/Params_'+fringe_field+'.dat', sep="\s+")

    particles_X = pd.read_table('Reference/X/Particles_'+fringe_field+'.dat', sep="\s+")
    particles_Y = pd.read_table('Reference/Y/Particles_'+fringe_field+'.dat', sep="\s+")
    particles_D = pd.read_table('Reference/D/Particles_'+fringe_field+'.dat', sep="\s+")

    tss_X = pd.read_table('Reference/X/TSS_analysis_'+fringe_field+'.dat', sep="\s+")
    tss_Y = pd.read_table('Reference/Y/TSS_analysis_'+fringe_field+'.dat', sep="\s+")
    tss_D = pd.read_table('Reference/D/TSS_analysis_'+fringe_field+'.dat', sep="\s+")

    temp_X = parameters_X.merge(particles_X, on=['I'])
    temp_Y = parameters_Y.merge(particles_Y, on=['I'])
    temp_D = parameters_D.merge(particles_D, on=['I'])
    data_X = temp_X.merge(tss_X, on=['I'])
    data_Y = temp_Y.merge(tss_Y, on=['I'])
    data_D = temp_D.merge(tss_D, on=['I'])

    #fig, ax = plot_3(data_X, 'CHROM_X', 'quadKx', 'SPIN-TUNE(1)', 'SPIN-TUNE(2)')
    #fig, ax = plot_3(data_Y, 'CHROM_Y', 'quadKy', 'SPIN-TUNE(1)', 'SPIN-TUNE(2)')
    #fig, ax = plot_3(data_D, 'ETA1', 'quadKz', 'SPIN-TUNE(1)', 'SPIN-TUNE(2)')

    #data_X.plot('CHROM_X', y=['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'])
    #data_Y.plot('CHROM_Y', y=['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'])
    #data_D.plot('ETA1', y=['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'])

    fig, ax = plot_3(data_X, 'CHROM_X', ['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'], 'SGx1', 'quadKx')
    fig, ax = plot_3(data_Y, 'CHROM_Y', ['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'], 'SGy1', 'quadKy')
    fig, ax = plot_3(data_D, 'ETA1', ['SPIN-TUNE(1)', 'SPIN-TUNE(2)','SPIN-TUNE(3)'], 'SGx1', 'quadKz')

plot_all('OFF')
plot_all('ON')
plt.show()
