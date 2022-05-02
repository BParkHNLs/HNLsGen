from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
(1.0,   1000.0, 2.30e-03, 50.94),
(1.0,    100.0, 7.90e-03, 7.71),
(1.0,     10.0, 9.49e-03, 6.49),
(1.5,   1000.0, 2.13e-03, 28.52),
(1.5,    100.0, 7.36e-03, 12.00),
(1.5,     10.0, 8.50e-03, 10.99),
(2.0,   1000.0, 2.05e-03, 43.36),
(2.0,    100.0, 7.12e-03, 12.12),
(2.0,     10.0, 7.99e-03, 10.76),
(3.0,   1000.0, 5.78e-03, 14.50),
(3.0,    100.0, 1.59e-02, 4.70),
(3.0,     10.0, 1.69e-02, 5.15),
(3.0,      1.0, 1.69e-02, 4.42),
(4.5,    100.0, 2.51e-02, 3.50),
(4.5,     10.0, 2.58e-02, 3.17),
(4.5,      1.0, 2.58e-02, 3.23),
(4.5,      0.1, 2.58e-02, 3.32),
]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=10000,muoneff=eff,displeff=1.0,timeevt=time,timejob=1,contingency=2)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







