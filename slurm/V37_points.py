from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
(1.0,   1000.0, 1.61e-03, 47.50),
(3.0,    100.0, 1.40e-02, 5.58),
(4.5,      1.0, 2.38e-02, 2.93),
(4.5,      0.1, 2.38e-02, 2.75),
]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=200000,muoneff=eff,displeff=1.0,timeevt=time,timejob=8,contingency=10)
  p.setConfig(cfg)
  points.append(p)


