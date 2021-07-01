from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [

(2.0,  10000.00, 1.87e-02, 18),
(2.0,   1000.00, 1.02e-01, 20),
(2.0,    100.00, 1.86e-01, 15),
(2.0,     10.00, 1.83e-01, 15),
(3.0,   1000.00, 8.48e-02, 15),
(3.0,    100.00, 1.57e-01, 16),
(3.0,     10.00, 1.59e-01, 22),
(3.0,      1.00, 1.59e-01, 22),
(4.5,    100.00, 3.60e-02, 19),
(4.5,     10.00, 3.64e-02, 20),
(4.5,      1.00, 3.64e-02, 19),
(4.5,      0.10, 3.64e-02, 19),
(5.5,     10.00, 1.31e-03, 48),
(5.5,      1.00, 1.31e-03, 51),
(5.5,      0.10, 1.38e-03, 50),
(5.5,      0.00, 1.35e-03, 43),


]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=5000,muoneff=eff,displeff=1.0,timeevt=time,timejob=12,contingency=3)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







