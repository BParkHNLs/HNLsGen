from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
#(3.0,   1000.0, 8.72e-02, 24.31),
(3.0,    100.0, 1.82e-01, 19.27),
(3.0,     10.0, 1.99e-01, 25.06),
#(3.0,      1.0, 1.86e-01, 19.98),
(4.5,    100.0, 1.30e-01, 21.08),
#(4.5,     10.0, 1.32e-01, 21.22),
(4.5,      1.0, 1.33e-01, 20.01),
#(4.5,      0.1, 1.31e-01, 20.66),
(5.5,     10.0, 1.39e-01, 24.12),
#(5.5,      1.0, 1.39e-01, 17.65),
(5.5,      0.1, 1.39e-01, 17.83),
#(5.5,      0.01, 1.40e-01, 17.71),

]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=150000,muoneff=eff,displeff=1.0,timeevt=time,timejob=6,contingency=6)
  p.setConfig(cfg)
  points.append(p)



