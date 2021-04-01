from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_vv_eff_time_s = [
(2.0,8.4e-05,2.0e-01,37),
(3.0,1.1e-05,1.5e-01,37),
(3.0,4.0e-04,1.7e-01,36),
(4.5,2.0e-04,3.5e-02,38),
#(4.5,6.0e-03,3.5e-02,39),
(5.5,1.0e-03,1.2e-03,75),
#(5.5,3.0e-02,1.2e-03,75),
]


points = []
for m,vv,eff,time in m_vv_eff_time_s:
  p   = Point(mass=m,ctau=None,vv=vv,ismaj=True)
  cfg = Config(nevtseff=100000,muoneff=eff,displeff=1.0,timeevt=time,timejob=35,contingency=3.)
  p.setConfig(cfg)
  points.append(p)

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







