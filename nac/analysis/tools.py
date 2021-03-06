__author__ = "Ivan Infante and Felipe Zapata"

# ================> Python Standard  and third-party <==========
import numpy as np
import os
import pyparsing as pa
from scipy.optimize import curve_fit

# ==================> Internal modules <==========
from nac.common import (hbar, r2meV, fs_to_cm)
# ==========================<>=================================


def autocorrelate(f):
    """
    Compute and returns the un-normalized and normalized autocorrelation
    of given function
    """
    d_f = f - f.mean()
    # Compute the autocorrelation function
    uacf = np.correlate(d_f, d_f, "full")[-d_f.size:] / d_f.size
    # Compute the normalized autocorrelation function
    nacf = uacf / uacf[0]
    return uacf, nacf


def gauss_function(x, sigma):
    """
    Gaussian function used to fit the dephasing time
    """
    return np.exp(-0.5 * (-x / sigma) ** 2)


def dephasing(f):
    """
    Computes the dephasing time of a given function using optical response
    formalisms:
    S. Mukamel, Principles of Nonlinear Optical Spectroscopy, 1995
    About the implementation we use the 2nd order cumulant expansion.
    See also eq. (2) in : Kilina et al. Phys. Rev. Lett., 110, 180404, (2013)
    To calculate the dephasing time tau we fit the dephasing function to a
    gaussian of the type : exp(-0.5 * (-x / tau) ** 2)
    """
    ts = np.arange(f.shape[0])
    cumu_ii = np.stack(np.sum(f[0:i]) for i in range(ts.size)) / hbar
    cumu_i = np.stack(np.sum(cumu_ii[0:i]) for i in range(ts.size)) / hbar
    deph = np.exp(-cumu_i)
    np.seterr(over='ignore')
    popt = curve_fit(gauss_function, ts, deph)[0]
    xs = np.exp(-0.5 * (-ts / popt[0]) ** 2)
    deph = np.column_stack((deph, xs))
    rate = popt[0]
    return deph, rate


def spectral_density(f):
    """
    Fourier Transform of a given function f using a dense grid with 100000 points.
    In the case of a FFT of a normalized autocorrelation function,
    this corresponds to a spectral density
    """
    f_fft = abs(1 / np.sqrt(2 * np.pi) * np.fft.fft(f, 100000)) ** 2
    # Fourier Transform of the time axis
    freq = np.fft.fftfreq(len(f_fft), 1)
    # Conversion of the x axis (given in cycles/fs) to cm-1
    freq = freq * fs_to_cm
    return f_fft, freq


def read_couplings(path_hams, ts):
    """
    This function reads the non adiabatic coupling vectors from the files
    generated for the NAMD simulations
    """
    files_im = [os.path.join(path_hams, 'Ham_{}_im'.format(i)) for i in range(ts)]
    xs = np.stack(np.loadtxt(fn) for fn in files_im)
    return xs * r2meV  # return energies in meV


def read_energies(path_hams, ts):
    """
    This function reads the molecular orbital energies of each state from
    the files generated for the NAMD simulations
    """
    files_re = [os.path.join(path_hams, 'Ham_{}_re'.format(i)) for i in range(ts)]
    xs = np.stack(np.diag(np.loadtxt(fn)) for fn in files_re)
    return xs * r2meV / 1000  # return energies in eV


def parse_list_of_lists(xs):
    """
    Parse a list of list of integers using pyparsing
    """
    enclosed = pa.Forward()  # Parser to be defined later
    natural = pa.Word(pa.nums)  # Natural Number
    # Nested Grammar
    nestedBrackets = pa.nestedExpr(pa.Suppress('['), pa.Suppress(']'), content=enclosed)
    enclosed << (natural | pa.Suppress(',') | nestedBrackets)
    try:
        rs = enclosed.parseString(xs).asList()[0]
        return list(map(lambda x: list(map(int, x)), rs))
    except pa.ParseException:
        raise RuntimeError("Invalid Macro states Specification")







