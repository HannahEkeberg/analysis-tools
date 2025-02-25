import curie as ci
import numpy as ni
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

class Count:

    def __init__(self):
        self.time_delay = [15*60,30*60, 60*60, 2*60*60, 3*60*60, 6*60*60, 10*60*60, 20*60*60, 30*60*60, 50*60*60, 100*60*60, 200*60*60, 500*60*60, 750*60*60, 1000*60*60]

    def efficiency(self, E):
        B0 = 5.131455901142433286e5
        B1 = 1.417943628682477097e1
        B2 = 4.844253332286091451e-2
        B3 = 1.218425547639006399e-5
        B4 = 2.519241334688014700e0
        return B0 * (np.exp(-B1*E**B2) ) * (1-np.exp(-B3 * E **B4))

    def t_count(self, lamb, I, t_delay, E):
        A0 = 3.7e4 # Bequerel ~ assuming micro curie
        Nc = 1e4 # Assuming desired counting time is 10 000
        countingTime = - np.log(  1 - (Nc*lamb / (A0 * I * np.exp(-lamb*t_delay) * self.efficiency(E)) ) ) / lamb
        logTerm = 1 - (Nc * lamb / (A0 * I * np.exp(-lamb * t_delay) * self.efficiency(E)))
        if logTerm != 1-5:
            countingTime = -np.log(logTerm) / lamb
        else:
            countingTime = np.nan
        return countingTime

    def decayConst(self, nucl):
        return ci.Isotope(nucl).decay_const()
    
    def halfLife(self, nucl):
        return ci.Isotope(nucl).half_life()

    def getCountingStatistics(self, isotopeDict):
        # isotopeDict --> {isotope: [[E1, I1], [E2, I2], ...]}
        countInfo = []
        for key, value in isotopeDict.items():
            isotope = key; E = value
            lamb = self.decayConst(isotope); thalf = self.halfLife(isotope)
            for t in self.time_delay:
                for e in range(len(E)):
                    countingTime = self.t_count(lamb=lamb, I=E[e][1], t_delay=t, E=E[e][0] )/60
                    countInfo.append([isotope, thalf/60, t/60, E[e][0], E[e][1]*100, countingTime])
        data = pd.DataFrame(
            countInfo, columns=['Isotope', 'Half life (min)', 'Time delay (min)', 'Energy', 'Intensity', 'Counting time (min)']).sort_values(by=['Half life (min)', 'Time delay (min)'])
        return data

    def saveCsv(self, dataframe, filename, directory=None):
        if directory==None:
            filename = os.getcwd() + '/' + filename + '.csv'
        else:
            filename = os.getcwd() + '/' + directory + '/' + filename + '.csv'
        print(filename)
        dataframe.to_csv(filename, index=False)

productsSc = {
    '63ZN': [[669.62, 8.2e-2], [962.06, 6.5e-2]],
    '45TI': [[719.6, 0.150e-2], [1408.1, 0.085e-2]],
    '43SC': [[372.9, 0.225]],
    '44SC': [[1157.9, 0.9989]],
    '62ZN': [[596.56, 0.26]],
    '42K': [[596.56, 0.1808]],
    '43K': [[372.76, 0.868]],
    'Sc-44m': [[271, 0.867]],
    '48V': [[983.525, 0.9998]],
    '44TI': [[78.32, 0.968]],
    '108AG': [[79.131, 0.066]],
    'Ag-105m': [[79.131, 0.066]],
    'Ag-108m': [[79.131, 0.066]]
    }
# dataSc = getCountingStatistics(productsSc)

productsAg = {
    '109CD': [[88.0336, 0.03644]],
    '107CD': [[93.124, 0.047], [828.943, 0.00163]],
    '104CD': [[66.6, 0.024], [83.5, 0.47], [559, 0.063], [625.7, 0.022], [709.2, 0.195]],
    '108AG': [[79.131, 0.066]],
    'Ag-108m': [[433.937, 0.905], [614.276, 0.898], [722.907, 0.908]],
    '106AG': [[873.46, 0.00199]],
    'Ag-106m': [[1045.83, 0.296], [717.34, 0.289]],
    '105AG': [[344.61, 0.42], [280.54, 0.31]],
    'Ag-105m': [[319.21, 0.0015]],
    '104AG': [[941.6, 0.25]],
    'Ag-104m': [[3034, 0.0023]],
    '103AG': [[118.74, 0.312], [148.20, 0.283]]
}

count = Count()

dataAg = count.getCountingStatistics(productsAg)
print(dataAg)
# saveCsv(dataframe=dataAg, filename='countingPlanSilver', directory=None)
# saveCsv(dataframe=dataSc, filename='countingPlanSilver', directory=None)