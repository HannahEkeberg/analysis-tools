import os
# from tools import *
import numpy as np
from scipy.interpolate import splev, splrep
import matplotlib.pyplot as plt

class Talys:

    def __init__(self, talysFilepath):
        self.talysFilepath = talysFilepath

    def generateTalysFiles(self, element, projectile, mass, energy, ldmodel=None, strength=None, astro=False, potential=None, outputfile=None):
        # ldmodel= 1,2,3,4...., strength=1,2,3,4,5
        if not os.path.exists(self.talysFilepath):
            os.makedirs(self.talysFilepath)

        inputfile = self.talysFilepath + '/inputfile.txt'
        inputfile_name = 'inputfile.txt'
        element = 'element ' + element + '\n'
        mass = 'mass ' + mass + '\n'
        projectile = 'projectile ' + projectile + '\n'
        energy = 'energy n0-' + energy + '.grid' + '\n'
        
        # if ldmodel==None:
        #     ldmodel = ''
        # else:    
        #     ldmodel = 'ldmodel ' + ldmodel + '\n'

        # if strength==None:
        #     strength = ''
        # else:    
        #     strength = 'strength ' + strength + '\n'

        # if astro==False:
        #     astro = ''
        # else:
        #     astro = 'astro y' + '\n'

        # if potential==None:
        #     potential = ''
        # else: 
        #     'jlmomp ' + potential[0] + '\n' + 'localomp' +  potential[1] + '\n'

        content = element + mass + projectile + energy #+ ldmodel + strength + astro
        with open( inputfile, 'w', encoding='utf-8') as file:
            file.write(content)
        
        os.chdir(self.talysFilepath)
        if outputfile == None:
            outputfile = 'outputfile.txt'
        command_string = 'talys <' + inputfile_name + '> ' + outputfile
        exit_code = os.system(command_string)
        os.chdir("..")

        if exit_code == 0:
            print("Talys ran successfully!" + self.talysFilepath)
        else:
            print(f"Talys failed, exit code: {exit_code}")

    def plotSpecificFile(self, label='TALYS-1.96', linestyle='-.', linewidth=1.5, color='orange', semilog_y=False):
        talysData= np.genfromtxt(self.talysFilepath)
        E = talysData[:,0]
        Cs = talysData[:,1]
        if semilog_y:
            plt.semilogy(E, Cs,label=label, color=color, linestyle=linestyle, linewidth=linewidth)
        else:
            plt.plot(E, Cs, label=label, linestyle=linestyle, color=color, linewidth=linewidth)
        return E,Cs

    def talysData(self, productZ, productA, targetFoil, isomerLevel = None):
        product = self.product(productZ, productA) #78, 198 --> 078193
        fileEnding = self.talysFileEnding(isomerLevel)
        filename = self.talysFilepath + targetFoil + '/rp' + product + fileEnding
        talysData = np.genfromtxt(filename)
        E = talysData[:,0]
        Cs = talysData[:,1]
        E, Cs = Tools_tendl().interpolate(E, Cs)
        return E, Cs

    def plotTalys(self,
    productZ,
    productA,
    targetFoil,
    isomerLevel = None,
    betaFeeding = None, # only beta+ beta-
    branchingRatio = None,
    parentIsomerLevel = None,
    ):
        E, Cs = self.talysData(productZ, productA, targetFoil, isomerLevel)
        if betaFeeding:
            E_parent, CsParent = self.correctForBetaFeeding(productZ, productA, targetFoil, betaFeeding, branchingRatio, parentIsomerLevel)
            Cs = Cs + CsParent
        plt.plot(E, Cs, label='TALYS-2.02', linestyle='-.', color='#ffe599')

    def correctForBetaFeeding(self, productZ, productA, targetFoil, betaFeeding, branchingRatio, parentIsomerLevel):
        if (betaFeeding  == 'beta+'):
            parentZ = str(int(productZ)+1); parentA = productA
        elif (betaFeeding == 'beta-'):
            parentZ = str(int(productZ)-1); parentA = productA
        E, Cs = self.talysData(parentZ, parentA, targetFoil, parentIsomerLevel)
        return E, Cs*branchingRatio

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

    def talysFileEnding(self, isomerLevel=None):
        return '.tot' if isomerLevel==None else '.L'+ isomerLevel


# talysFilePath = os.getcwd() + '/../talys_v2.04/'
# Talys(talysFilePath).plotTalys('78', '193', 'Ir', '05')


class Tools_tendl:

    def interpolate(self, x, y, Elimit=None, zeroPadding=False):
        if Elimit==None:
            Elimit=200
        if zeroPadding:
            x, y = self.zeroPadding(x,y)
        tck = splrep(x, y, s=0)
        x_new = np.linspace(1, Elimit, 1000)
        y_new = splev(x_new, tck, der=0)
        return x_new, y_new

    def zeroPadding(self, x, y):
        if x[0]!=0:
            zero_padding = np.linspace(0,x[0]-0.5,10)
            zeros_y = np.zeros((len(zero_padding)))
            x = np.concatenate((zero_padding, x))
            y = np.concatenate((zeros_y, y))
        return x, y