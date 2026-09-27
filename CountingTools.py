import curie as ci
import numpy as ni
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from scipy.constants import e, N_A
from nuclearanalysistools.Tendl import Tendl
from nuclearanalysistools.pandastools import PandasTools

class Count:

    def __init__(self, foil, target, numberOfTargetNuclei, beamEnergy, beamParticle):
        self.foil = foil # 'SC'
        self.target = target # {'Sc45': 1.0}
        self.numberOfTargetNuclei = numberOfTargetNuclei # ask chat gpt for estimate, eg. cu = 2.12e24 nuclei/m2
        self.beamEnergy = beamEnergy #'70'
        self.beamParticle = beamParticle # 'proton
        self.time_delay = [0, 15*60,30*60, 60*60, 2*60*60, 3*60*60, 6*60*60, 10*60*60, 20*60*60, 30*60*60, 50*60*60, 100*60*60, 200*60*60, 500*60*60, 750*60*60, 1000*60*60]
    
    def efficiency(self, E):
        # B0 = 5.131455901142433286e5
        # B1 = 1.417943628682477097e1
        # B2 = 4.844253332286091451e-2
        # B3 = 1.218425547639006399e-5
        # B4 = 2.519241334688014700e0
        B0 = 1.776545259072640664e+45
        B1 = 1.034688117827637086e+02
        B2 = 7.750520399491661432e-03
        B3 = 4.648389369558355438e-03
        B4 = 1.241889791865289361e+00
        return B0 * (np.exp(-B1*E**B2) ) * (1-np.exp(-B3 * E **B4))

    def A0(self, product, productZ, productA, isomerLevel=None):#, areal_density, thickness_foil):    #63ZN, '30, '63'
        # Assuming thickness foil = 25 micro m
        molar_mass = ci.Element(self.foil).mass # amu
        beam_current = 100 * e #nA * elementary charge --> protons/s
        decay_const = self.decayConst(product)
        t_irr = 2*3600 #(s)
        E, Cs = Tendl(self.target, self.beamParticle).tendlData(productZ=productZ, productA=productA)#, isomerLevel=None, Elimit = None)
        indeces = np.where( (E > int(self.beamEnergy) - 1) & (E < int(self.beamEnergy) + 1))
        Cs_av = np.mean(Cs[indeces])*1e-3 #barn
        E_av = np.mean(E[indeces])
        A0 = Cs_av*beam_current*self.numberOfTargetNuclei*(1-np.exp(-decay_const*t_irr)) #[Bq]
        return A0

    def estimatedA(self, A0, decay_const, t_delay):
        A = A0 * np.exp(-decay_const * t_delay) # e^{-\lambda t_delay}. If we want to check counting stats after eob. 
        A = f"{A:.2e}"
        return A

    def N_T(self, areal_density, molar_mass):   ### Not currently in use.....
        # areal density --> g/m2, molar_mass --> g/mol, avogadros number (N_T) --> 1/mol
        return N_A * (areal_density/molar_mass) # nuclei/m^2


    def t_count(self, A0, lamb, I, t_delay, E):
        A0 = 3.7e4 # Bequerel ~ assuming micro curie
        Nc = 1e4 # Assuming desired counting time is 10 000
        countingTime = - np.log(  1 - (Nc*lamb / (A0 * I * np.exp(-lamb*t_delay) * self.efficiency(E)) ) ) / lamb
        return countingTime

    def decayConst(self, isotope):
        return ci.Isotope(isotope).decay_const()
    
    def halfLife(self, isotope):
        return ci.Isotope(isotope).half_life()

    def getCountingStatistics(self, isotopeDict):
        # isotopeDict --> {isotope: [[[E1, I1], [E2, I2], ...], Z, A, isomerLevel]}
        countInfo = []
        for key, value in isotopeDict.items():
            isotope = key; isotopeValues = value
            E = isotopeValues[0]; isotopeZ = isotopeValues[1]; isotopeA=isotopeValues[2]
            lamb = self.decayConst(isotope); thalf = self.halfLife(isotope)
            A0 = self.A0(isotope, isotopeZ, isotopeA)
            for t in self.time_delay:
                for e in range(len(E)):
                    countingTime = self.t_count(A0, lamb=lamb, I=E[e][1], t_delay=t, E=E[e][0] )/60
                    A = self.estimatedA(A0, lamb, t)
                    countInfo.append([isotope, thalf/60, t/60, A, self.foil, E[e][0], E[e][1]*100, countingTime])
        data = pd.DataFrame(
            countInfo, columns=['Isotope', 'Half life (min)', 'Time delay (min)', 'Estimated A (Bq)', 'Foil', 'Energy', 'Intensity', 'Counting time (min)']).sort_values(by=['Half life (min)', 'Time delay (min)'])
        return data
    