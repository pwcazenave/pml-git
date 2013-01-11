""" Quick test of doing the FFT of a time series to get the harmonic amplitude and phase """

import numpy as np
import matplotlib.pyplot as plt
import datetime

from scipy.stats import mode
from matplotlib.dates import date2num


def plot_spectrum(t, Fs):
    """
    Plot a spectrum of time series t.

    Lifted from
    http://glowingpython.blogspot.co.uk/2011/08/how-to-plot-frequency-spectrum-with.html

    Parameters
    ----------

    t : ndarray
        Time series of surface elevations.
    Fs : float
        Sampling rate (number of samples per time).

    """

    n = len(t) # length of the signal
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range

    Y = np.fft.fft(t)/n # fft computing and normalization
    Y = Y[range(n/2)]

    plt.plot(frq, abs(Y), 'r') # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')


tt = np.genfromtxt('/users/modellers/pica/Data/NTSLF/formatted/AVO.txt', delimiter=',')

t = tt[645000:745000, :-2] # some relatively clean time series
tbar = np.mean(t[~np.isnan(t[:, -1]), -1])
t[np.isnan(t[:, -1])] = tbar # replace NaNs with mean
t[:, -1] = t[:, -1] - tbar # remove the mean
tz = t[:, -1] # tide elevations

# Create a time series in decimal
# TODO: Use MJD.
td = []
for i, dt in enumerate(t[:, :6]):
    dtt = [int(x) for x in dt]
    temp_t = date2num(datetime.datetime(dtt[0], dtt[1], dtt[2], dtt[3], dtt[4], dtt[5]))
    if temp_t < 70000:
        print temp_t, dtt, i
    td.append(temp_t)

td = np.asarray(td)

ts = mode(np.diff(td))[0] # use the mode the diffs to get the sampling rate (per day)
fs = 1.0/ts

plt.subplot(2,1,1)
plt.plot(td, tz)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.subplot(2,1,2)
plot_spectrum(tz, fs)
plt.show()
