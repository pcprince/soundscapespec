# soundscapespec
Python implementation of R function which returns a kHz binned spectrum as described by Kasten et al. (2012) for the description of a soundscape.

The soundscape frequency spectrum is based on the computation of a spectrogram power spectral density using Welch'smethod (Welch & June, 1967). Based on the R implementation included as [part of the seewave package](http://rug.mnhn.fr/seewave/HTML/MAN/soundscapespec.html).

Argument       | Description
---------------|-------------
`wave`         | Audio samples
`f`            | Sample rate
`wl`           | FFT window length
`normalised`   | Whether or not result is normalised between 0 and 1

Returns a two-column numeric matrix, the first column returning the frequency (kHz) bands and the second column returning the power value within each frequency band.

## References

Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote environmental assessment laboratory's acoustic library: an archive for studying soundscape ecology. Ecological Informatics, 12, 50-67.

Welch, P.D., June (1967). The use of the fast Fourier transform for the estimation of power spectra: a method based on time-averaging over short, modified periodograms. IEEE Transactions on Audio and Electroacoustics, 15: 70-73. 
