from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
(3.0,184.,3.93e-03,149),
]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=200000,muoneff=eff,displeff=1.0,timeevt=time,timejob=20,contingency=3.)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







