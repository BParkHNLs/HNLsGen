from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float


m_vv_eff_time_s = [
(0.5,8.6e-02,6.0e-03,104),
(1.0,2.7e-03,3.3e-03,160),
(1.5,3.5e-04,2.0e-03,216),
(2.0,8.4e-05,5.7e-04,435),
(3.0,1.1e-05,4.7e-03,149),
#(3.0,4.0e-04,5.6e-03,107),
(4.5,2.0e-04,4.8e-04,700),
#(4.5,6.0e-03,4.8e-04,700),
]

points = []
for m,vv,eff,time in m_vv_eff_time_s:
  p   = Point(mass=m,ctau=None,vv=vv,ismaj=True)
  cfg = Config(nevtseff=300000,muoneff=eff,displeff=1.0,timeevt=time,timejob=35,contingency=3.)
  p.setConfig(cfg)
  points.append(p)

  # 

  #print('mass={:.1f} vv={:.1e} before={:.1e}, after={:.1e}'.format(p.mass,p.vv,displEff[m][vv],getCtauEff(p.ctau*4.,1000.))) # 4 is the assumed beta*gamma factor...
  #cfg.stamp()
  







