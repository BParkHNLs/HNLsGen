import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate

mycolors = ['blue', 'darkcyan', 'lawngreen','orange', 'red']

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])


# old filter eff settings, A-M is using these
m_ctau_eff_time_s = [

(1.0,  10000.0, 2.24e-04, 367),
(1.0,   1000.0, 1.79e-03, 68),
(1.0,    100.0, 5.72e-03, 36),
(1.0,     10.0, 6.57e-03, 38),
(1.5,  10000.0, 1.07e-04, 732),
(1.5,   1000.0, 8.83e-04, 107),
(1.5,    100.0, 3.09e-03, 40),
(1.5,     10.0, 3.56e-03, 42),
(2.0,  10000.0, 4.40e-05, 1663),
(2.0,   1000.0, 3.66e-04, 229),
(2.0,    100.0, 1.34e-03, 79),
(2.0,     10.0, 1.57e-03, 67),
(3.0,   1000.0, 1.93e-03, 62),
(3.0,    100.0, 5.15e-03, 39),
(3.0,     10.0, 5.46e-03, 41),
(3.0,      1.0, 5.47e-03, 30),
(4.5,    100.0, 3.53e-04, 215),
(4.5,     10.0, 4.59e-04, 163),
(4.5,      1.0, 4.59e-04, 167),
(4.5,      0.1, 4.58e-04, 196),

]

# super new ones, new filter efficiency settings
#m_ctau_eff_time_s = [
#
#(1.0,  10000.00, 1.69e-04, 427),
#(1.0,   1000.00, 1.34e-03, 70),
#(1.0,    100.00, 4.24e-03, 37),
#(1.0,     10.00, 4.86e-03, 29),
#(1.5,  10000.00, 1.02e-04, 678),
#(1.5,   1000.00, 8.42e-04, 98),
#(1.5,    100.00, 2.71e-03, 46),
#(1.5,     10.00, 2.99e-03, 49),
#(2.0,  10000.00, 4.68e-05, 1565),
#(2.0,   1000.00, 3.67e-04, 196),
#(2.0,    100.00, 1.25e-03, 75),
#(2.0,     10.00, 1.42e-03, 65),
#(3.0,   1000.00, 2.04e-03, 52),
#(3.0,    100.00, 4.82e-03, 33),
#(3.0,     10.00, 5.02e-03, 30),
#(3.0,      1.00, 5.02e-03, 30),
#(4.5,    100.00, 3.66e-04, 207),
#(4.5,     10.00, 4.39e-04, 165),
#(4.5,      1.00, 4.42e-04, 143),
#(4.5,      0.10, 4.42e-04, 131),
#
#]

masses = [1.0, 1.5, 2.0, 3.0, 4.5]

ctaus = {}
effs  = {}
fits  = {}
popts = {}
tcks  = {}

fout = open('effs_V31.py', 'w')
fout.write('effs={}\n')


for mass in masses:
  ctaus[mass] = []
  effs[mass] = []
  #fits[mass] = [] only one element
  for m,ctau,eff,time in reversed(m_ctau_eff_time_s):
    if m==mass:
      ctaus[mass].append(ctau*1E-03)
      effs[mass].append(eff)  


# spline fit degree 1 (~piecewise interpolation between two points)
for mass in masses:
  #fit = np.polyfit(ctaus[mass],effs[mass],5)
  #fit_fn = np.poly1d(fit)
  #print fit_fn(0.1)
  #fits[mass] = fit_fn
  #popt,pcov = optimize.curve_fit(piecewise_linear, ctaus[mass],effs[mass])
  #popts[mass] = popt
  #print np.array(ctaus[mass])
  #print np.array(effs[mass])
  tck = interpolate.splrep(np.array(ctaus[mass]),np.array(effs[mass]), k=1, s=0)
  tcks[mass] = tck


# plot the results and print them to file
fig, ax0 = plt.subplots()
ax0.set_xlabel('c$\\tau$ (m)')
ax0.set_ylabel('filter Eff')
#ax0.set_ylim(1e-06,1.)
ax0.set_yscale('log')
ax0.set_xscale('log')
ax0.grid(which='both', axis='both')
for i,mass in enumerate(masses):
  ctau_max = max(ctaus[mass])
  ctau_min = min(ctaus[mass])
  big_ctau = np.logspace(start=np.log10(ctau_min), stop=np.log10(ctau_max), num=100)
  ax0.plot(ctaus[mass],effs[mass], 'o', label='mass={:.1f}GeV'.format(mass), color=mycolors[i])
  #fit_fn = fits[mass]
  #ax0.plot(big_ctau,fit_fn(big_ctau), label='fit mass={:.1f}GeV'.format(mass), color=mycolors[i])
  #p = popts[mass]
  #ax0.plot(big_ctau, piecewise_linear(big_ctau, *p), color=mycolors[i])
  tck = tcks[mass]
  ax0.plot(big_ctau, interpolate.splev(big_ctau, tck, der=0), color=mycolors[i])

  fout.write( 'effs[{m}] = {{}}\n'.format(m=mass))
  fout.write( 'effs[{m}]["ctau"] = {a}\n'.format(m=mass, a=list(big_ctau)))
  fout.write( 'effs[{m}]["eff"] = {a}\n'.format(m=mass, a=list(interpolate.splev(big_ctau, tck, der=0))))


ax0.legend(loc='upper right', fontsize='x-small', numpoints=1)
fig.savefig('effs.pdf')
fout.close()

  




