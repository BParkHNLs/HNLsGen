from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
(1.0,   1000.0, 2.30e-03, 50.94),
(1.0,     10.0, 9.49e-03, 6.49),
(1.5,    100.0, 7.36e-03, 12.00),
(1.5,     10.0, 8.50e-03, 10.99),
(2.0,    100.0, 7.12e-03, 12.12),
(2.0,     10.0, 7.99e-03, 10.76),
(3.0,    100.0, 1.59e-02, 4.70),
(3.0,     10.0, 1.69e-02, 5.15),
(4.5,    100.0, 2.51e-02, 3.50),
(4.5,      1.0, 2.58e-02, 3.23),
]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=400000,muoneff=eff,displeff=1.0,timeevt=time,timejob=6,contingency=5)
  p.setConfig(cfg)
  points.append(p)


