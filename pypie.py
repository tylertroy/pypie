import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import re
from matplotlib import patches
from itertools import cycle
from matplotlib.markers import MarkerStyle
import matplotlib.ticker as plticker
from os.path import join, split
from os import listdir
from pylab import connect
from pylab import searchsorted
from math import sqrt
from numpy import average
from numpy import mean
from numpy import linspace
from numpy import vstack
from numpy import std
from collections import deque

""" Function Definitions """

def load(filepath=''):
    with open(filepath, 'r') as f:
        energy = [float(elem) for elem in f.readline().split()[1:]]
        current = [float(elem) for elem in f.readline().split()[1:]]
        data = list(zip(*[[int(elem) for elem in column] 
            for column in [line.split() for line in f.readlines()]]))
        time, counts = data[0], data[1:]
    return energy, current, time, counts
def cal_mass(m1, m2, t1, t2):
    """ 
    Calculate the proportionality, k, and time zero, t0, constants
    given two masses and their time/data points for t = k*sqrt(mass)
    """
    # m1, m2 = masses
    # t1, t2 = times
    t0 = -1 * ((-t2 + t1) / (sqrt(m1) - sqrt(m2))) * sqrt(m1) + t1
    k = (-t2 + t1) / (sqrt(m1) - sqrt(m2))
    return k, t0
def map_mass(k, t0, npoints):
    """ """
    return [index_to_mass(t, k, t0) for t in xrange(npoints)]
def index_to_mass(index, k, t0):
    """ Returns mass at a given point """
    return ((index-t0)/k)**2
def mass_to_index(mass, k, t0):
    """ Returns mass at a given point """
    return int(k*sqrt(mass) + t0)
def save_mass(data_dict):
    filepath = ''.join((data_dict['path'], data_dict['name'][4:], '_MZ')) 
    counts2d = [['Energy'] + data_dict['energy']] + \
               [['Current'] + data_dict['current']] + \
                 zip(*([data_dict['mz']] + data_dict['counts']))
    data_string = '\n'.join(('\t'.join((str(item) for item in row)) for row in counts2d))
    with open(filepath, 'w') as f:
        f.write(data_string)
def avg_pie(center=1, width=2, counts=[[]], axis=0):
    """
    Return the average of each point between center +/_ width
    """
    halfwidth = int(round(width * 0.5))
    start, stop = (center - halfwidth, center + halfwidth)
    return average(a=counts[start:stop], axis=axis)
def find_peaks(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    from numpy import NaN, Inf, arange, isscalar, asarray, array

    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return maxtab

""" Class Definitions """
        
class Pie(object):
    """A 3D data processing library for mass spectrometric data series.
    """
    _marker = cycle(MarkerStyle.markers.keys())
    def __init__(self, filepath):
        super(Pie, self).__init__()
        self.load(filepath)
        self.mz = {}
        self.pie = {}
        self.xax = []
    def pie_colors(self, space, limits):
        self.colors = iter(space(linspace(*limits)))
    def get_peaks(self, index, mincount):
        self.peak_energy = index
        index = index if type(index) is int else self.energy.index(index)
        self.peaks = find_peaks(self.counts[index], mincount, self.mass)
    def peaks_as_ratio(self):
        counts = list(zip(*self.peaks))[1]
        total_counts = sum(counts)
        self.peaks = [ (mass, count / total_counts ) for mass, count in self.peaks ]
    def save_peaks(self, filename):
        header = ''.join(['m/z', '\t','counts at ', str(self.peak_energy), 'eV','\n']) 
        with open(filename, 'w') as f:
            f.write(header)
            for peak, count in self.peaks:
                string = ''.join([str(peak), '\t', str(count), '\n'])
                f.write(string)
    def load(self, filepath):
        self.filepath, self.filename = split(filepath)
        self.energy, self.current, self.time, self.counts = load(filepath)
        self.energy = [round(e, 3) for e in self.energy]
        self.xbins = len(self.energy)
        self.ybins = len(self.time)
        self.xpts = range(self.xbins)
        self.ypts = range(self.ybins)
        self.counts_sum = np.sum(np.array(self.counts).T, axis=1)
    def stdev(self, files=()):
        pass
    def save(self, filepath):
        if filepath:
            path, name = split(filepath)
            if not path:
                path = self.filepath
        filepath = join(path, name)
    def by_mass(self, key):
        """ Sort function to sort by PIE mass center.
        """
        return self.pie[key]['mz']
    def ms_cursor(self, index=None, block=True):
        if index:
            try:
                index = index if type(index) is int else self.energy.index(index)
                energy = index
            except ValueError:
                print('{} is not an available energy\nTry one of these:'.format(index))
                print(self.energy)
                return
            y = self.counts[index]
        else:
            y = self.counts_sum
            energy = '{}-{}'.format(self.energy[0], self.energy[-1])
        x = self.ypts       
        fig, ax = plt.subplots()
        c = SnaptoCursor(ax, x, y)
        connect('motion_notify_event', c.mouse_move)
        ax.plot(x, y, 'r')
        ax.set_xlabel('Point')
        ax.set_ylabel('Ion Count')
        ax.set_title('Mass Spectrum at {}eV'.format(energy))
        plt.xlim(0,max(x))
        plt.show(block=block)
    def ms_calibrate(self, m1, m2, t1, t2):
        terms = m1, m2, t1, t2
        self.k, self.t0 = cal_mass(*terms)
        self.mass = [index_to_mass(t, self.k, self.t0) 
            for t in range(self.ybins)]
    def ms_plot(self, index=None):
        if index:
            index = index if type(index) is int else self.energy.index(index)
            energy = index
            spectrum = self.counts[index]
        else:
            spectrum = self.counts_sum
            energy = '{}-{}'.format(self.energy[0], self.energy[-1])
        plt.plot(self.mass, spectrum)
        plt.xlim((min(self.mass), max(self.mass)))
        plt.xlabel('m/z', fontsize=20)
        plt.ylabel('Ion Counts', fontsize=20)
        plt.title('Mass Spectrum at {}eV'.format(energy))
        plt.show()
    def ms_save(self, index=None, path=''):
        if index:
            index = index if type(index) is int else self.energy.index(index)
            counts = self.counts[index]
            energy = self.energy[index]
        else:
            counts = self.counts_sum
            energy = '{}-{}'.format(self.energy[0], self.energy[-1])
        if not path:
            energy_string = str(energy).replace('.','p')
            path = '{}_{}_MS.txt'.format(self.filename, energy_string)

        header = ''.join(['m/z', '\t','counts at ', str(energy), 'eV','\n']) 
        with open(path, 'w') as f:
            f.write(header)
            for mass, count in zip(self.mass, counts):
                string = '\t'.join([str(mass), str(count)]) + '\n'
                f.write(string)
    def ms_save_all(self, filename=''):
        header = 'm/z\t' + '\t'.join([ str(energy) for energy in self.energy ]) + '\n'
        with open(filename, 'w') as f:
            f.write(header)
            for mass, count in zip(self.mass, zip(*self.counts)):
                counts =  '\t'.join([str(c) for c in count])
                string = ''.join([str(mass), '\t', counts, '\n'])
                f.write(string)
    def pie_slice(self, center, width, label=None):
        label = label if label else center
        mcenter = center
        halfwidth = width * 0.5
        mstart, mstop = center - halfwidth, center + halfwidth
        (center, start, stop) = (mass_to_index(mcenter, self.k, self.t0),
                                 mass_to_index(mstart, self.k, self.t0),
                                 mass_to_index(mstop, self.k, self.t0))
        counts = list(zip(*self.counts))
        avg = mean(counts[start:stop], axis=0)
        pie_slice = np.sum(counts[start:stop], axis=0)
        pie = { 'slice'   : pie_slice,
                'ptcoord' : (center, start, stop),
                'mzcoord' : (mcenter, mstart, mstop),
                'ptslice' : (start, abs(start - stop), center),
                'mzslice' : (mstart, abs(mstart - mstop), mcenter),
                'marker'  : next(self._marker),
                'mz'      : mcenter
              }
        self.pie[label] = pie
    def pie_save(self, pie_keys=None, path=''):
        if not path:
            path = '{}_PIEs'.format(self.filename)
        if not pie_keys:
            pie_keys = list(self.pie.keys())
            try:
                pie_keys = sorted(pie_keys, key=self.by_mass)
                pie_keys = [ 'mz'+key]
            except:
                pass
        if all([ key in self.pie.keys() for key in pie_keys ]):
            data = []
            for key in pie_keys:
                data.append(self.pie[key]['slice'])
            data.insert(0, self.energy)
            data = list(zip(*data))
            data.insert(0, ['Energy /ev'] + pie_keys)
            savestring = '\n'.join([ '\t'.join([str(element) 
                for element in line]) 
                    for line in data ])
            with open(path, 'w') as f:
                f.write(savestring)
        else:
            for key in pie_keys:
                if key not in self.pie.keys():
                    print('%s is not a valid PIE key. show available keys with object.pie.keys()' % key)
                    return
    def pie_normalize(self, pie_keys=None):
        if not pie_keys:
            pie_keys = list(self.pie.keys())
        for key in pie_keys:
            pie = self.pie[key]['slice']
            self.pie[key]['slice'] = [ count-min(pie) for count in pie ]
            pie = self.pie[key]['slice']
            self.pie[key]['slice'] = [ count/max(pie) for count in pie ]
    def pie_current_correction(self, pie_keys=None):
        if not pie_keys:
            pie_keys = list(self.pie.keys())
        current = [ c/max(self.current) for c in self.current ]
        for key in pie_keys:
            pie = self.pie[key]['slice']
            self.pie[key]['slice'] = [ (count / c) for count, c in zip(pie, current) ]
    def pie_show_slices(self, index=None, filepath=None, xlim=None, ylim=None, params=None):
        """ Visualize the data slices extracted by self.pie_slice. If no
        filepath is supplied, spectrum will not be saved. """
        if index:
            index = index if type(index) is int else self.energy.index(index)
            spectrum = self.counts[index]
            energy = self.energy[index]
        else:
            spectrum = self.counts_sum
            energy = '{}-{}'.format(self.energy[0], self.energy[-1])
        default = {
            'color':'k',
        }
        params = params if params else {}
        params.update(default)
        slices = [ pie['mzslice'] for pie in self.pie.values() ]
        fig, ax = plt.subplots()
        fig.set_size_inches(18.5,10.5)
        loc = plticker.MultipleLocator(base=10.0)
        plt.minorticks_on()
        ax.xaxis.set_major_locator(loc)
        ax.plot(self.mass, spectrum, **params)
        for start, width, center in slices:
            ymax = ax.viewLim.ymax
            r = patches.Rectangle((start, 0), width, ymax, color='b', fill=True, alpha=0.2)
            ax.add_artist(r)
            ax.vlines(center, 0, ymax, color='r')
        xlim = xlim if xlim else (min(self.mass), max(self.mass))
        ylim = ylim if ylim else (min(spectrum), max(spectrum))
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('m/z', fontsize=16)
        plt.ylabel('Ion Counts', fontsize=16)
        plt.title('PIE Slices Over Mass Spectrum at {}eV'.format(energy))
        if filepath:
            path, name = split(filepath)
            if not path:
                path = self.filepath
            plt.savefig(filepath)
        else:
            plt.show()
    def pie_delete(self, key):
        """ a simple wrapper for the self.pie dict pop function """
        self.pie.pop(key)
    def pie_plot(self, pie_keys=None, params=None):
        params = params if params else {}
        numbers = re.compile(r'(\d+)')
        if not pie_keys:
            pie_keys = sorted(list(self.pie.keys()), key=self.by_mass)
        nspec = len(pie_keys)
        self.pie_colors(matplotlib.cm.hot, (0,0.7,nspec))
        fig, ax = plt.subplots()
        fig.set_size_inches(18.5,10.5)
        for key, color in zip(pie_keys, self.colors):
            label = numbers.sub(r'$_{\1}$', key)
            marker = self.pie[key]['marker']
            pie = self.pie[key]['slice']
            ax.plot(self.energy, pie, color=color, marker=marker, label=label, markersize=8)
        xlim =  (min(self.energy), max(self.energy))
        plt.xlim(xlim)
        plt.xlabel('Energy (eV)', fontsize=16)
        plt.ylabel('Intensity (counts)', fontsize=16)
        plt.title('Photoionization Efficiency Curve')
        plt.legend(loc=2, borderaxespad=0., fontsize='x-large')
        plt.show()
    def pie_background_correction(self, background_key, pie_keys=None):
        """ Subtract PIE labelled "background_key" from all other PIES.
        """
        if not pie_keys:
            pie_keys = list(self.pie.keys())
        bg_key = background_key
        bg = self.pie[bg_key]['slice']
        for key in pie_keys:
            pie = self.pie[key]['slice']
            self.pie[key]['slice'] = [ count - b for count, b in zip(pie, bg) ]
    def pie_normalize_to(self, index, pie_keys=None):
        """ Normalize extracted PIEs to the value specified at index.
        """
        index = index if type(index) is int else self.energy.index(index)
        if not pie_keys:
            try: 
                pie_keys = sorted(list(self.pie.keys()), key=self.by_mass)
            except:
                pie_keys = sorted(list(self.pie.keys()))
        for key in pie_keys:
            pie = self.pie[key]['slice']
            normal = pie[index]
            self.pie[key]['slice'] = [ count / normal for count in pie ]
    def pie_info(self):
        info = ['slice', 'marker']
        for key in self.pie.keys():
            print('{}:'.format(key))
            for key, val in self.pie[key].items():
                if key not in info:
                    print('{}: {}'.format(key, val))
            print()
    def current_plot(self):
        plt.plot(self.energy, self.current)
        plt.xlabel('Energy (eV)')
        plt.ylabel('Current (A)')
        plt.title('Beam Current Measured by Kiethley-6485')
        plt.show()

class SnaptoCursor:
    """
    Like Cursor but the crosshair snaps to the nearest x,y point
    For simplicity, I'm assuming x is sorted
    """
    def __init__(self, ax, x, y):
        self.ax = ax
        self.lx = ax.axhline(color='k')  # the horiz line
        self.ly = ax.axvline(color='k')  # the vert line
        self.x = x
        self.y = y
        # text location in axes coords
        self.txt = ax.text( 0.7, 0.9, '', transform=ax.transAxes)

    def mouse_move(self, event):

        if not event.inaxes: return

        x, y = event.xdata, event.ydata

        indx = searchsorted(self.x, [x])[0]
        x = self.x[indx]
        y = self.y[indx]
        # update the line positions
        self.lx.set_ydata(y )
        self.ly.set_xdata(x )

        self.txt.set_text( 'x=%1.2f, y=%1.2f'%(x,y) )
        # print ('x=%1.2f, y=%1.2f'%(x,y))
        plt.draw()

if __name__ == '__main__':
    pass
