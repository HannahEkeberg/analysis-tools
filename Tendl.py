# from tools import *
import numpy as np
import requests
import matplotlib.pyplot as plt
from urllib.request import urlopen
from scipy.interpolate import splev, splrep

class Tendl_Tools:

    def interpolate(self, x, y, Emin, Elimit, zeroPadding=False):
        if Elimit==None:
            Elimit=200
        if Emin == None:
            Emin=1
        if zeroPadding:
            x, y = self.zeroPadding(x,y)
        tck = splrep(x, y, s=0)
        x_new = np.linspace(Emin, Elimit, 1000)
        y_new = splev(x_new, tck, der=0)
        
        return x_new, y_new

    def zeroPadding(self, x, y):
        if x[0]!=0:
            zero_padding = np.linspace(0,x[0]-0.5,10)
            zeros_y = np.zeros((len(zero_padding)))
            x = np.concatenate((zero_padding, x))
            y = np.concatenate((zeros_y, y))
        return x, y


class Tendl:

    def __init__(self, target, beamParticle):
        self.target = target # target = {"Ir191": 0.373, "Ir193": 0.627}, {"Ta181": 0.9998799, "Ta180m": 0.0001201}
        self.beamParticle = beamParticle
        self.color = '#a6376e'
        # self.color = '#b61163'
        self.linestyle = '--'
        self.label='TENDL-2025'
        self.linewidth=1.0


    def get_q_value(self, productZ, productA, isomerLevel):
        print(self.target)
        targetFoil = list(self.target.keys())[0][0:2]
        product = self.product(productZ, productA)
        fileEnding = self.tendlFileEnding(isomerLevel)
        print("Z: " + productZ + ", A: " + productA)
        for t in self.target.keys():
            q_val = self.retrieveTendlQvalue(self.tendlUrl(targetFoil, t, product, fileEnding))
            if not isinstance(q_val, str):
                q_val = '%.1f' %q_val
            print(t + ': ' + q_val)
        pass

    def tendlDeuteronData(self, productZ, productA, isomerLevel = None):
        targetFoil = list(self.target.keys())[0][0:2]
        product = self.product(productZ, productA)
        fileEnding = self.tendlFileEnding(isomerLevel)
        E = []
        Cs = []
        for t in self.target.keys():
            data = self.retrieveTendlDataFromUrl(
                self.tendDeuteronlUrl(targetFoil, t, product, fileEnding), t
            )
            if isinstance(data[0], np.ndarray):
                E.append(data[0])
                Cs.append(data[1])
        if len(E)==0 or len(Cs)==0:
            "TENDL: No data found for target: " + targetFoil + " for productZ" + productZ + "and product A: " + productA
            raise Exception
        CsSummed = sum(Cs)
        E = E[0]
        E, Cs = Tendl_Tools().interpolate(E, CsSummed)
        return E, Cs

    def tendlData(self, productZ, productA, isomerLevel=None,Emin = None, Elimit=None):
        # targetFoil = list(self.target.keys())[0][0:1] # Just setting 1 for Y
        print(type(self.target))
        target = list(self.target.keys())[0]
        if target[1].isdigit():
            targetFoil = target[0:1]
        else:
            targetFoil = target[0:2]
        product = self.product(productZ, productA)
        fileEnding = self.tendlFileEnding(isomerLevel)
        E = []
        Cs = []
        for t in self.target.keys():
            data = self.retrieveTendlDataFromUrl(
                self.tendlUrl(targetFoil, t, product, fileEnding), t
            )
            if isinstance(data[0], np.ndarray):
                E.append(data[0])
                Cs.append(data[1])
        if len(E)==0 or len(Cs)==0:
            "TENDL: No data found for target: " + targetFoil + " for productZ" + productZ + "and product A: " + productA
            # raise Exception
        try:
            CsSummed = sum(Cs)
            E = E[0]
        except:
            min_len = min(len(c) for c in Cs)
            CsSummed = sum(c[:min_len] for c in Cs)
            E = E[0][:min_len]
        E, Cs = Tendl_Tools().interpolate(E, CsSummed,Emin, Elimit)
        return E, Cs

    def plotTendl23(self, productZ, productA, isomerLevel = None): #, feeding = None, branchingRatio = None, parentIsomerLevel = None):
        # try:
        E, Cs = self.tendlDeuteronData(productZ, productA, isomerLevel)
        # if feeding == 'beta+' or feeding == 'beta-':
            # CsParent = self.correctForFeeding(productZ, productA, feeding, branchingRatio, parentIsomerLevel)[1]
            # Cs = Cs + CsParent
        plt.plot(E, Cs, label='TENDL-2025', linestyle='--', color='royalblue')
    # except:
        # print("Unable to retrive tendl data, perhaps no internet connection?")

    def plotTendl23Unique(self, productZ, productA, isomerLevel = None, color=None, lineStyle=None, label=None, semilog_y=False):
        # try:
        if color==None:
            color='#741b47'
        if lineStyle==None:
            lineStyle='--'
        if label==None:
            label = 'TENDL-2025'
        E, Cs = self.tendlData(productZ, productA, isomerLevel)
        if semilog_y:
            plt.semilogy(E, Cs,label=label, color=color, linestyle=linestyle, linewidth=linewidth)
        else:
            plt.plot(E, Cs, label=label, linestyle=lineStyle, color=color, linewidth=1)
        # except:
        # print("Unable to retrive tendl data, perhaps no internet connection?")

    def plotTendl23Unique_feeding(self, productZ, productA, isomerLevel = None,  
                                  betaPlusDecayChain = None, betaMinusDecayChain = None, isomerDecayChain = None, 
                                  color=None, lineStyle=None, label=None, semilog_y=False):
        # {isotope: [productZ, branchingRatio isomerLevel]} #beta+/beta-
        # {isotope: [branchingRatio isomerLevel]} #isomer
        try:
            if color==None:
                color = self.color
            if lineStyle==None:
                lineStyle = self.linestyle
            if label==None:
                label = self.label
            E, Cs = self.tendlData(productZ, productA, isomerLevel)
            Cs_betaplus = []; Cs_betaminus = []; Cs_isomer = []
            if betaPlusDecayChain:
                for i in list(betaPlusDecayChain.keys()):
                    Z = betaPlusDecayChain[i][0]
                    branchingRatio= betaPlusDecayChain[i][1]
                    isomerLevel = betaPlusDecayChain[i][2]
                    E_bp, Cs_bp = self.tendlData(Z, productA, isomerLevel)
                    Cs_betaplus.append(Cs_bp*branchingRatio)
            if betaMinusDecayChain:
                for i in list(betaMinusDecayChain.keys()):
                    Z = betaMinusDecayChain[i][0]
                    branchingRatio= betaMinusDecayChain[i][1]
                    isomerLevel = betaMinusDecayChain[i][2]
                    E_bm, Cs_bm = self.tendlData(Z, productA, isomerLevel)
                    Cs_betaminus.append(Cs_bm*branchingRatio)
            if isomerDecayChain:
                for i in list(isomerDecayChain.keys()):
                    branchingRatio= isomerDecayChain[i][0]
                    isomerLevel = isomerDecayChain[i][1]
                    E_i, Cs_i = self.tendlData(productZ, productA, isomerLevel)
                    Cs_isomer.append(Cs_i*branchingRatio)
            totCs = Cs + sum(Cs_betaplus) + sum(Cs_betaminus) + sum(Cs_isomer)
            # plt.plot(E, totCs, label='TENDL-2023', linestyle='--', color='royalblue')
            plt.plot(E, totCs, label=label, linestyle=lineStyle, linewidth=self.linewidth,  color=color)
        except:
            print("Unable to retrive tendl data, perhaps no internet connection?")

    def plotdataWithMultipleFeeding(self, productZ, productA, isomerLevel, betaPlusDecayChain = None, betaMinusDecayChain = None, isomerDecayChain = None):
        # {isotope: [productZ, branchingRatio isomerLevel]} #beta+/beta-
        # {isotope: [branchingRatio isomerLevel]} #isomer
        try:
            E, Cs = self.tendlDeuteronData(productZ, productA, isomerLevel)
            Cs_betaplus = []; Cs_betaminus = []; Cs_isomer = []
            if betaPlusDecayChain:
                for i in list(betaPlusDecayChain.keys()):
                    Z = betaPlusDecayChain[i][0]
                    branchingRatio= betaPlusDecayChain[i][1]
                    isomerLevel = betaPlusDecayChain[i][2]
                    E_bp, Cs_bp = self.tendlDeuteronData(Z, productA, isomerLevel)
                    Cs_betaplus.append(Cs_bp*branchingRatio)
            if betaMinusDecayChain:
                for i in list(betaMinusDecayChain.keys()):
                    Z = betaMinusDecayChain[i][0]
                    branchingRatio= betaMinusDecayChain[i][1]
                    isomerLevel = betaMinusDecayChain[i][2]
                    E_bm, Cs_bm = self.tendlDeuteronData(Z, productA, isomerLevel)
                    Cs_betaminus.append(Cs_bm*branchingRatio)
            if isomerDecayChain:
                for i in list(isomerDecayChain.keys()):
                    branchingRatio= isomerDecayChain[i][0]
                    isomerLevel = isomerDecayChain[i][1]
                    E_i, Cs_i = self.tendlDeuteronData(productZ, productA, isomerLevel)
                    Cs_isomer.append(Cs_i*branchingRatio)
            totCs = Cs + sum(Cs_betaplus) + sum(Cs_betaminus) + sum(Cs_isomer)
            plt.plot(E, totCs, label='TENDL-2023', linestyle='--', color='royalblue')
        except:
            print("Unable to retrive tendl data, perhaps no internet connection?")

    def product(self, productZ, productA):
        if len(productZ) <= 2:
            productZ = '0' + productZ
        else:
            productZ = productZ
        if len(productA) <= 2:
            productA = '0' + productA
        else:
            productA = productA
        return productZ + productA

    def tendDeuteronlUrl(self, targetFoil, target, product, fileEnding):
        if len(target)<5:
            target = target[0:2] + '0' + target[2:]
        return ('https://tendl.web.psi.ch/tendl_2023/deuteron_file/'
        + targetFoil + '/' + target
        + '/tables/residual/rp'
        + product + fileEnding)

    def tendlUrl(self, targetFoil, target, product, fileEnding):
        if self.beamParticle == 'deuteron':
            beam_file = 'deuteron_file/'
        elif self.beamParticle == 'proton':
            beam_file = 'proton_files/'
        elif self.beamParticle == 'alpha':
            beam_file = 'alpha_file/'
        elif self.beamParticle == 'neutron':
            beam_file = 'neutron_files/'
        elif self.beamParticle == 'He3':
            beam_file = 'he3_file/'
        else:
            raise Exception("Invalid beam particle. Must be deuteron or proton. Was: " + self.beamParticle)
        target = self.formatTargetLength(targetFoil, target)
        # url= ('https://tendl.web.psi.ch/tendl_2023/' 
        url= ('https://tendl.imperial.ac.uk/tendl_2025/' 
            #   https://tendl.imperial.ac.uk/tendl_2025/proton_files/Ta/Ta181/tables/residual/rp073177.tot
            #   https://tendl.imperial.ac.uk/tendl_2025/proton_file/Ta/Ta181/tables/residual/rp073177.tot
            + beam_file
            + targetFoil + '/' + target
            + '/tables/residual/rp'
            + product + fileEnding)
        print(url)
        return url
        # return ('https://tendl.web.psi.ch/tendl_2023/' 
        # + beam_file
        # + targetFoil + '/' + target
        # + '/tables/residual/rp'
        # + product + fileEnding)

    def formatTargetLength(self, targetFoil, targetIsotope):
        # Cu65 --> Cu065. Ir193=Ir193
        isotopeNumber = targetIsotope[len(targetFoil):]
        number = ''.join(filter(str.isdigit, isotopeNumber))
        isomer = isotopeNumber[len(number):] 
        formattedIsotopeNumber = number if len(number)==3 else '0' + number
        return targetFoil + formattedIsotopeNumber + isomer

    def tendlFileEnding(self, isomerLevel=None):
        return '.tot' if isomerLevel==None else '.L' + isomerLevel

    def retrieveTendlDataFromUrl(self, url, target):
        try:
            tendlData = requests.get(url).text.split("\n")[27:] # skipping 27 first lines in tendl file
            tendlData = np.genfromtxt(tendlData)
            abundance = self.target[target]
            E = tendlData[:,0]
            Cs = tendlData[:,1]
            return E, Cs*abundance
        except:
            print('Unable to retrieve tendlData from url: ' + url)
            return 0,0
    
    def retrieveTendlQvalue(self, url):
        try:
            tendlData = requests.get(url).text.split("\n")
            Q_value = float([l for l in tendlData if "Q-value" in l][0].split(":")[1].split()[0])
            return Q_value
        except:
            print('Unable to retrieve Q value and data from url: ' + url)
            return "Undefined"

    # def retrieveDataFromUrlWithNumpy(self, url):
    #     tendl_data = np.genfromtxt(urlopen("https://tendl.web.psi.ch/tendl_2023/deuteron_file/Ir/Ir193/tables/residual/rp078193.L05"), delimiter=" ")
    #     energy = tendl_data[:,0]
    #     xs = tendl_data[:,1]
    #     return energy, xs

    # def mapFeedingObject(self, feeding):
    #     new = []
    #     for item in feeding:
    #         new.append(vars(item))
    #     # print(new)
            



