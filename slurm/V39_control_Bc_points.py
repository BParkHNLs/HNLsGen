from python.common  import Point,Config,getCtauEff

# List of points to be generated
#   - N.B.: mass must be a float

p = Point(mass=999,ctau=999,vv=999)
cfg = Config(nevtseff=100000,muoneff=2.072e-01,displeff=1.0,timeevt=100,timejob=6,contingency=3.)
p.setConfig(cfg)

points = []
points.append(p)








