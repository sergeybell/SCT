import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#
#    textstr = '\n'.join((
#    r'$\mu=%.2f$' % (mu, ),
#    r'$\mathrm{median}=%.2f$' % (median, ),
#    r'$\sigma=%.2f$' % (sigma, )))

def plot_3_rays(data, params, X, Y1, Y2, Y3, num_rays):
    print(params['I'])
    textstr = '\n'.join((
    r'$\mathrm{CASE}=%.2f$'                    % (params['I'], ),
    r'$\xi_x=%.2f$'             % (params['CHROM_X'], ),
    r'$\xi_y=%.2f$'             % (params['CHROM_Y'], ),
    r'$\eta_1=%.2f$'             % (params['ETA1'], )
    ))

    fig, ax = plt.subplots(3,1,sharex=True)
    t = data.loc[data['ray'] == 1, ['iteration']]
    for i in range (num_rays):
        ax[0].plot(t, data.loc[data['ray'] == i+1, [Y1]]); ax[0].set_ylabel(Y1);
        ax[1].plot(t, data.loc[data['ray'] == i+1, [Y2]]); ax[1].set_ylabel(Y2);
        ax[2].plot(t, data.loc[data['ray'] == i+1, [Y3]]); ax[2].set_ylabel(Y3);
    ax[2].set_xlabel(X)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax[0].text(0.05, 0.9, textstr, transform=ax[0].transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
    return fig, ax

fringe_field = 'OFF'

CASE = 5

parameters_X = pd.read_table('X/Params_'+fringe_field+'.dat', sep="\s+")
#parameters_Y = pd.read_table('Y/Params_'+fringe_field+'.dat', sep="\s+")
#parameters_D = pd.read_table('D/Params_'+fringe_field+'.dat', sep="\s+")
parameters_X = parameters_X.loc[parameters_X['I'] == CASE]
#parameters_Y = parameters_Y.loc[parameters_Y['I'] == CASE]
#parameters_D = parameters_D.loc[parameters_D['I'] == CASE]

#print(parameters_X)

CASE = str(CASE)
colnames_orbital=['iteration', 'ray', 'X', 'A', 'Y', 'B', 'T', 'D']
TRPRAY_X = pd.read_table('X/TRPRAY_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_orbital)
#TRPRAY_Y = pd.read_table('Y/TRPRAY_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_orbital)
#TRPRAY_D = pd.read_table('D/TRPRAY_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_orbital)

colnames_spin=['iteration', 'ray', 'S_X', 'S_Y', 'S_Z']
TRPSPI_X = pd.read_table('X/TRPSPI_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_spin)
#TRPSPI_Y = pd.read_table('Y/TRPSPI_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_spin)
#TRPSPI_D = pd.read_table('D/TRPSPI_'+fringe_field+'_'+CASE+'.dat', sep="\s+", comment='#', names = colnames_spin)

#print(TRPRAY_X)
#print(TRPSPI_X)

tracking_data_X= TRPRAY_X.merge(TRPSPI_X, on=["iteration","ray"])
#tracking_data_Y= TRPRAY_Y.merge(TRPSPI_Y, on=["iteration","ray"])
#tracking_data_D= TRPRAY_D.merge(TRPSPI_D, on=["iteration","ray"])

#print(tracking_data_X)

fig, ax = plot_3_rays(tracking_data_X, parameters_X, 'iteration', 'S_X', 'S_Y', 'S_Z', 4)
#fig, ax = plot_3_rays(tracking_data_Y, parameters_Y, 'iteration', 'S_X', 'S_Y', 'S_Z', 4)
#fig, ax = plot_3_rays(tracking_data_D, parameters_D, 'iteration', 'S_X', 'S_Y', 'S_Z', 4)

plt.show()
