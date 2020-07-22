from scipy.signal import get_window
from scipy.fft import fft
from scipy.io import wavfile

import numpy as np
import math

def stdft(wave, f, wl, steps, half_wl):

    W = get_window("hamming", wl)

    z = []

    for i in range(len(steps)):
        s = steps[i]
        f = fft(wave[s:s+wl]*W)
        # Keep only the relevant frequencies (half of the FFT)
        z.append(f[0:half_wl])

    z = np.transpose(z)

    # Scale by the original number of fft points
    z = z / wl
    # Multiply by 2 to save the total energy (see http://www.ni.com/white-paper/4278/en/, section Converting from a Two-Sided Power Spectrum to a Single-Sided Power Spectrum)
    z = 2 * abs(z)

    return z

def soundscapespec(wave, f, wl, normalised):

    half_wl = int(wl/2)

    # Welch spectrogram

    n = len(wave)

    # Positions of the sliding window for the STFT
    steps = range(0, n-wl, half_wl)
    # Number of records processed
    n_recs = len(steps)

    # Frequency values of the spectrum
    x = np.linspace(f/wl, (f/2) - (f/wl), half_wl)/1000

    # STFT with the parameters set above
    y = stdft(wave, f, wl, steps, half_wl)
    # Square of the summed spectrum
    y = pow(y, 2)
    # Sum of successive spectra
    y = np.sum(y, axis=1)
    # Normalize for number of records and frequency
    y = y / (n_recs * f)

    # Frequency bins

    # Truncate frequency digits
    freq = np.trunc(x)

    # Eliminate frequency values between 0 and 1 kHz
    freq = freq[freq!=0]
    # Eliminate amplitude values between 0 and 1 kHz
    spec = y[np.trunc(x)!=0]

    # 1 kHz frequency bins
    bins = np.unique(freq)

    bin_lengths = [0] * len(bins)

    for i in bins:

        index = int(i) - 1
        bin_lengths[index] = len(freq[freq==i])

    # Width of each frequency bin
    bin_width = 1000 * wl / f

    # Remove small frequency bin at the extreme right of the spectrum
    bin_out = []

    for bi in range(len(bin_lengths)):

        if bin_lengths[bi] < np.trunc(bin_width):

            bin_out.append(bi)

    if len(bin_out) > 0:

        bins = np.delete(bins, bin_out)

    # Compute the energy of each frequency bin
    val = [0] * len(bins)

    for i in bins:

        index = int(i) - 1
        val[index] = sum(spec[freq == bins[index]])

    if normalised:
        # Compute the vector length
        vlen = math.sqrt(sum(np.square(val)))
        # Normalize to the range [0,1]
        val = np.array(val) / vlen

    return np.column_stack((bins, val))
