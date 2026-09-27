from scipy.interpolate import splev, splrep
import numpy as np
from datetime import datetime

class Tools:

    def interpolate(self, x, y, xlimit=None, zeroPadding=False):
        if xlimit == None:
            xlimit = 40
        if zeroPadding:
            x, y = self.zeroPadding(x,y)
        tck = splrep(x, y, s=0)
        x_new = np.linspace(1, xlimit, 1000)
        y_new = splev(x_new, tck, der=0)
        return x_new, y_new

    def zeroPadding(self, x, y):
        if x[0]!=0:
            zero_padding = np.linspace(0,x[0]-0.5,10)
            zeros_y = np.zeros((len(zero_padding)))
            x = np.concatenate((zero_padding, x))
            y = np.concatenate((zeros_y, y))
        return x, y

    def colors(self):
        # return ['mediumpurple', 'cyan', 'palevioletred', 'darkorange', 'forestgreen', 'orchid', 'dodgerblue', 'lime', 'crimson', 'indianred']
        # minty #B4C9C4, 
        # pistachio #93C572
        return [ 'crimson', 'dodgerblue', 'forestgreen', 'palevioletred', 'sienna', 'indianred','darkgoldenrod', 'firebrick', 'orchid',  'darkorange','dodgerblue', 'lime','mediumpurple', '#a4c483', 'pink', '#B4C9C4', '#93C572', 'cyan',]
    
    def date_diff(self, start_time, end_time, units=None):
        #t1 = '09/23/2025 20:19:12';t2 = '09/25/2025 15:45:00'
        fmt = '%m/%d/%Y %H:%M:%S'

        dt1 = datetime.strptime(start_time, fmt)
        dt2 = datetime.strptime(end_time, fmt)
        delta = dt2 - dt1
        # delta = dt1 - dt2
        if units==None:
            return delta  # days h s
        if units=='d':
            return delta.days
        elif units=='h':
            return delta.hours
        elif units=='min':
            return delta.min
        else:
            raise Exception('Invalid unit: ' + units)
        
