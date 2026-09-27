from scipy.interpolate import splev, splrep
import numpy as np
import os

class PandasTools:

    def concatDataframes(self, dataframes=[]):
        return pd.concat(datafranes, ignore_index=True)
        

    def saveCsv(self, dataframe, filename, directory=None):
        if directory==None:
            filename = os.getcwd() + '/' + filename + '.csv'
        else:
            filename = os.getcwd() + '/' + directory + '/' + filename + '.csv'
        print(filename)
        dataframe.to_csv(filename, index=False)