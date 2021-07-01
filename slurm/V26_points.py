from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_ctau_eff_time_s = [
#(1.0,184.  ,3.97e-03,50),
#(1.0,1000. ,2.40e-03,100),
(1.0,10000.,2.40e-04,1000),
]

points = []
for m,ctau,eff,time in m_ctau_eff_time_s:
  p   = Point(mass=m,ctau=ctau,vv=None,ismaj=True)
  cfg = Config(nevtseff=20000,muoneff=eff,displeff=1.0,timeevt=time,timejob=12,contingency=3.)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







