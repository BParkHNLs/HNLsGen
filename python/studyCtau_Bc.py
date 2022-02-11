'''
Plotting script for the 3D ctau plot as a function of VV and mass - Bc grid
Includes ctau levels and a grid of points that were requested for central production
Figure 8 of AN
'''

from common import getCtau

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import LogLocator

#plt.style.use('./mystyle.mplstyle')

#masses = np.array( [float('{:.2f}'.format(0.5+i*0.1)) for i in range(0,56) ] )# GeV
masses = np.array( [float('{:.2f}'.format(0.8+i*0.1)) for i in range(0,53) ] )# GeV
# only have a few numbers for QCD_corrections 

#masses = np.logspace( np.log10(0.5), np.log10(6.0),  56, base=10)
#vvs = np.logspace( np.log10(1e-06), np.log10(3e-04), 50, base=10)
vvs = np.logspace( np.log10(1e-06), np.log10(1), 200, base=10)
#vvs = np.logspace( np.log10(1e-06), np.log10(1), 10, base=10)
#ctau_levels = np.array([0.00001,0.0001,0.001, 0.01, 0.1, 1.,10.])
ctau_levels = np.array([0.00001,0.0001,0.001, 0.01, 0.1, 1.])
mycolors = ['black','blue', 'darkcyan', 'green', 'lawngreen','yellow','orange', 'red']

ctaus = []
for vv in vvs:
  ctaus_m = []
  for m in masses:
    ctaus_m.append(getCtau(m,vv)/1000)
  ctaus.append(ctaus_m)

ctaus = np.array(ctaus)
masses, vvs = np.meshgrid(masses, vvs)

fig, ax0 = plt.subplots()
pcm = ax0.pcolor(masses, vvs, ctaus,
                   norm=colors.LogNorm(vmin=ctaus.min(), vmax=ctaus.max()),
                   cmap='PiYG_r')
cbar = fig.colorbar(pcm, ax=ax0, extend='max')
cbar.set_label('c$\\tau$ (m)')
CS = ax0.contour(masses, vvs, ctaus, levels=ctau_levels, colors=mycolors) #cmap='PiYG_r', norm=colors.LogNorm(vmin=ctaus.min(), vmax=ctaus.max()))
ax0.set_xlabel('mass (GeV)')
ax0.set_ylabel('$V^2$')
ax0.set_ylim(1e-06,1.)
ax0.set_xlim(0.8,6)
ax0.set_yscale('log')
ax0.set_xscale('log')
ax0.grid(which='both', axis='both')


xs = {}
ys = {}

for i,ctau in enumerate(ctau_levels):
  print ctau
  xs[ctau] = []
  ys[ctau] = []

  p = CS.collections[i].get_paths()[0]
  v = p.vertices
  for j,m in enumerate([2.0,3.,4.5,5.5]):
    #print m
    res = np.where(v == m)
    #print res
    if res and res[0].size != 0:
      mass = v[ res[0][0],0]
      vv   = v[ res[0][0],1]
      if vv < 0.01 and getCtau(mass,vv)<=(1000+500) and mass > 2. : # only consider points below 10^{-2} and 
        xs[ctau].append( mass )
        ys[ctau].append( vv )

  ax0.plot(xs[ctau],ys[ctau], 'o', label='ctau={:.1e}m'.format(ctau), color=mycolors[i])
ax0.legend(loc='upper right', fontsize='x-small', numpoints=1)

fig.savefig('ctau_Bc.pdf')
fig.savefig('ctau_Bc.png')



