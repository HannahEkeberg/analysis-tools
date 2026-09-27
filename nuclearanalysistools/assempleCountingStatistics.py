from CountingTools import *




count_Ge = Count(foil='GE', target={'Ga69': 0.60108, 'Ga71': '0.39892'}, numberOfTargetNuclei=2.56e20, beamEnergy=55, beamParticle='proton')

count_Ge.getCountingStatistics()