import numpy as np
import obspy
from obspy import read
import scipy.interpolate
import scipy.signal
import const
import sys
import io


class Sac(object):
    def __init__(self, file_name):
        try:    
            self.sac = read(file_name)[0]
            stats = self.sac.stats.sac
        except:
            f = open("errors1.txt", "a")
            f.write('sac_read_error : %s\n' %file_name)
            f.close()
            raise RuntimeError('{} : cannot read sac file'.format(file_name))

        nt = 86400 * self.sac.stats.sampling_rate
        if stats.npts < nt:
            f = open("errors1.txt", "a")
            f.write('short_data_length : %s\n' %file_name)
            f.close()
            raise RuntimeError('{} : Data points is not enough, {}'.format(file_name, stats.npts))
        if np.all(self.sac.data == 0):
            f = open("errors1.txt", "a")
            f.write('empty_sac_file : %s\n' %file_name)
            f.close()
            raise RuntimeError('{} : No data'.format(file_name))

        self.sac.filter('bandpass', freqmin=const.freq_min, freqmax=const.freq_max, corners=const.corner, zerophase=True)
        self.sac.decimate(factor=5, no_filter=True)
        self.sac.data = self.sac.data**2
        self.sac.filter('lowpass', freq=const.low_pass, corners=const.corner, zerophase=True)
        self.sac.decimate(factor=20, no_filter=True)

        self.x = [stats.stla, stats.stlo, stats.stel / 1000.0]
        self.path = file_name


class Splines(object):
    def __init__(self, v):
        x = np.arange(v.shape[0])
        y = np.arange(-v.shape[1]//2+1, v.shape[1]//2+1)
        f = scipy.interpolate.RectBivariateSpline(x, y, v)
        self.a = np.zeros(v.shape + (4,), dtype=np.float32)
        self.a[:, :, 0] = f(x, y)
        self.a[:, :, 1] = f(x, y, dy=1)
        self.a[:, :, 2] = f(x, y, dy=2) / 2
        self.a[:, :, 3] = f(x, y + 1) - np.sum(self.a, axis=2)


    def __call__(self, x, y):
        y1, y2 = np.modf(y + self.a.shape[1]//2)
        y2 = y2.astype(np.int32)
        v = self.a[x, y2, 0].copy()
        for i in range(1, 4):
            v += pow(y1, i) * self.a[x, y2, i]
        return v


    def grad(self, x, y):
        y1, y2 = np.modf(y + self.a.shape[1]//2)
        y2 = y2.astype(np.int32)
        v = self.a[x, y2, 1].copy()
        v += 2 * self.a[x, y2, 2] * y1
        v += 3 * self.a[x, y2, 3] * y1 * y1
        return v


class CrossCorrelations(object):
    def __init__(self, v, i_s, j_s, i_c, j_c, i_p, j_p, file_format):
        v = np.array(v, dtype=np.float32)
        self.i_s = np.array(i_s, dtype=np.int32)
        self.j_s = np.array(j_s, dtype=np.int32)
        self.i_c = np.array(i_c, dtype=np.int32)
        self.j_c = np.array(j_c, dtype=np.int32)
        self.i_p = i_p
        self.j_p = j_p
        
        if v.shape[0] < const.min_pair_num:
            return
        
        v = v[:, v.shape[1]//2-40:v.shape[1]//2+40+1]
        self._f = Splines(v)
        self.removed = np.zeros(v.shape[0], dtype=np.bool_)


    def __call__(self, x, y):
        return self._f(x, y)


    def grad(self, x, y):
        return self._f.grad(x, y)