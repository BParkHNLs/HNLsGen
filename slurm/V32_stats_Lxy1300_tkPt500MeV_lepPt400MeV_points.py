from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [

(1.0,  10000.00, 1.69e-04, 427),
(1.0,   1000.00, 1.34e-03, 70),
(1.0,    100.00, 4.24e-03, 37),
(1.0,     10.00, 4.86e-03, 29),
(1.5,  10000.00, 1.02e-04, 678),
(1.5,   1000.00, 8.42e-04, 98),
(1.5,    100.00, 2.71e-03, 46),
(1.5,     10.00, 2.99e-03, 49),
(2.0,  10000.00, 4.68e-05, 1565),
(2.0,   1000.00, 3.67e-04, 196),
(2.0,    100.00, 1.25e-03, 75),
(2.0,     10.00, 1.42e-03, 65),
(3.0,   1000.00, 2.04e-03, 52),
(3.0,    100.00, 4.82e-03, 33),
(3.0,     10.00, 5.02e-03, 30),
(3.0,      1.00, 5.02e-03, 30),
(4.5,    100.00, 3.66e-04, 207),
(4.5,     10.00, 4.39e-04, 165),
(4.5,      1.00, 4.42e-04, 143),
(4.5,      0.10, 4.42e-04, 131),

]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=5000,muoneff=eff,displeff=1.0,timeevt=time,timejob=12,contingency=3.)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







