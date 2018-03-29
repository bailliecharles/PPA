read('partitions.nex')   # read nexus partition file
read('AA.phy')           # read in alignment (phylip)
a = var.alignments[0]    # assign a the alignment
a.setCharPartition('PartitionFinder') # and the partitions

d = Data()               # attach Data
read('*con.tre')         # read in MrBayes consensus tree
t=var.trees[0]           # assign to t
t.data=d                 # assign the data to the tree

### Define partitions. MUST be the same as the MrBayes run.
### Here that is seven partitions (zero-based) as defined by PartitionFinder2

t.newComp(partNum=0, free=0, spec='empirical')
t.newRMatrix(partNum=0, free=0, spec='mtREV24')
t.setPInvar(partNum=0, free=1, val=0.1)
t.setNGammaCat(partNum=0, nGammaCat=4)
t.newGdasrv(partNum=0, free=1, val=0.5)

t.newComp(partNum=1, free=0, spec='empirical')
t.newRMatrix(partNum=1, free=0, spec='mtREV24')
t.setPInvar(partNum=1, free=0, val=0)
t.setNGammaCat(partNum=1, nGammaCat=4)
t.newGdasrv(partNum=1, free=1, val=0.5)

t.newComp(partNum=2, free=0, spec='empirical')
t.newRMatrix(partNum=2, free=0, spec='mtREV24')
t.setPInvar(partNum=2, free=1, val=0.1)
t.setNGammaCat(partNum=2, nGammaCat=4)
t.newGdasrv(partNum=2, free=1, val=0.5)

t.newComp(partNum=3, free=0, spec='empirical')
t.newRMatrix(partNum=3, free=0, spec='mtREV24')
t.setPInvar(partNum=3, free=0, val=0.1)
t.setNGammaCat(partNum=3, nGammaCat=4)
t.newGdasrv(partNum=3, free=1, val=0.5)

t.newComp(partNum=4, free=0, spec='empirical')
t.newRMatrix(partNum=4, free=0, spec='mtREV24')
t.setPInvar(partNum=4, free=1, val=0.1)
t.setNGammaCat(partNum=4, nGammaCat=4)
t.newGdasrv(partNum=4, free=1, val=0.5)

t.newComp(partNum=5, free=0, spec='empirical')
t.newRMatrix(partNum=5, free=0, spec='mtREV24')
t.setPInvar(partNum=5, free=1, val=0.1)
t.setNGammaCat(partNum=5, nGammaCat=4)
t.newGdasrv(partNum=5, free=1, val=0.5)

t.newComp(partNum=6, free=0, spec='empirical')
t.newRMatrix(partNum=6, free=0, spec='mtREV24')
t.setPInvar(partNum=6, free=1, val=0.1)
t.setNGammaCat(partNum=6, nGammaCat=4)
t.newGdasrv(partNum=6, free=1, val=0.5)

# This bit is necessary, can't remember why
t.data = d
t.calcLogLike()

# slurp up the posterior samples 
ps = PosteriorSamples(t, runNum=2, program='mrbayes', mbBaseName='AA.nex', verbose=2)

# make sure you import random 
import random

# Start loop. Here a thousand draws are randomly made from the post burn-in (20%)
# set of parameter estimates

for sampNum in random.sample(range(2000, 10000), 1000):
     t2 = ps.getSample(sampNum)
     t2.data = d
     t2.simulate()
     t2.data.alignments[0].writePhylip(fName='mySims.phy', interleave=False, flat=True, append=True)

