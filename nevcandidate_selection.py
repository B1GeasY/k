#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 10:59:56 2020

@author: k
"""

###this code is for selecting [NeV] emission objects, compute SNR and simply fit a gaussian profile to find the linecenter

from astropy.io import fits
from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt
from specutils import Spectrum1D, SpectralRegion
import astropy.units as u
from astropy.nddata import StdDevUncertainty
from specutils.fitting import fit_lines
from astropy.modeling import models
from specutils.fitting import fit_generic_continuum
from specutils.manipulation import extract_region
from specutils.manipulation import noise_region_uncertainty
from specutils.manipulation import box_smooth
from specutils.analysis import gaussian_fwhm,centroid
import os
from specutils.analysis import snr

from kapteyn import kmpfit
from astropy.table import Table


def onegauss(xval, pp):
    """The single Gaussian model used to fit the emission lines 
    Parameter: the scale factor, central wavelength in logwave, line FWHM in logwave
    """
    
    term1 = np.exp( - (xval - pp[1])**2 / (2. * pp[2]**2) )
    yval = pp[0] * term1 / (np.sqrt(2.*np.pi) * pp[2])
    return yval

def _residuals(pp,data):
    xval, yval, weight = data
    return (yval - onegauss(xval,pp))/weight






def fit_line(name,wave,flux,err,z):
    
    rest_wave = wave
    uncertainty = StdDevUncertainty(err*u.Jy)
    spec1 = Spectrum1D(spectral_axis = rest_wave*u.AA,flux=flux*u.Jy,uncertainty=uncertainty)
    spec = box_smooth(spec1, width=3)
    nev_region = SpectralRegion(3395*u.AA,3408*u.AA)+SpectralRegion(3440*u.AA,3456*u.AA)
    
    nev_lineregion = SpectralRegion(3415*u.AA,3435*u.AA)
    
    nev_spec = extract_region(spec,SpectralRegion(3390*u.AA,3460*u.AA))
    
    line_spec = extract_region(nev_spec,nev_lineregion)
#    nev_spec = extract_region(spec,nev_region)
#    SNR1 = snr(line_spec)   ###median snr
    
    line_cont = fit_generic_continuum(nev_spec,exclude_regions = nev_lineregion)
    cont_fit = line_cont(nev_spec.spectral_axis)
    
#    peak1 = centroid(nev_spec) 
    norm_nev_spec = Spectrum1D(spectral_axis = nev_spec.spectral_axis,flux = nev_spec.flux/cont_fit)
    sub_nev_spec = Spectrum1D(spectral_axis = nev_spec.spectral_axis,flux = nev_spec.flux-cont_fit)
    
    line_fit = kmpfit.Fitter(residuals = _residuals,data = (sub_nev_spec.spectral_axis/u.AA,\
                                                            sub_nev_spec.flux/u.Jy,\
                                                                nev_spec.uncertainty.array))
    
    line_fit_ini = [0,3426,3]
    line_fit.parinfo = [{'limits':(0.,10.**10)},{'limits':(3416,3436)},{'limits':(0,6)}]
    line_fit.fit(params0=line_fit_ini)
    snrregion = SpectralRegion((line_fit.params[1]-2*line_fit.params[2])*u.AA,(line_fit.params[1]+2*line_fit.params[2])*u.AA)
#    print(snrregion)
    calc_spectrum = extract_region(nev_spec, snrregion)
    
    flux = calc_spectrum.flux
    uncertainty = calc_spectrum.uncertainty.array * nev_spec.uncertainty.unit
    SNR = np.median(flux / uncertainty, axis=-1)
    
    peak = line_fit.params[1]
    
    wave = nev_spec.spectral_axis/u.AA
#    print(snr(nev_spec,snrregion).value)
    plt.figure()
    plt.plot(nev_spec.spectral_axis,nev_spec.flux)
    plt.plot(nev_spec.spectral_axis,cont_fit,'orange')
    plt.plot(nev_spec.spectral_axis,nev_spec.uncertainty.array,'grey')
    plt.plot(sub_nev_spec.spectral_axis,sub_nev_spec.flux,'green')
    plt.axvline(3426.)
    plt.plot(nev_spec.spectral_axis/u.AA,onegauss(nev_spec.spectral_axis/u.AA,line_fit.params),'r')
    plt.savefig('/Users/k/Documents/MMT_spec/MMT_LockmanHole/NeV/24b_NeV/'+name+'.eps',format='eps',dpi=300)
#    plt.close()
    
    return SNR,line_fit.params[1]
    
    
    

fits_path = '/Users/k/Documents/MMT_spec/MMT_LockmanHole/24b_spec/'
fitslist = os.listdir(fits_path)
fitslist.remove('.DS_Store')

SNR = []
Peak = []
name = []

for i in range(len(fitslist)):
    hdu = fits.open(fits_path+fitslist[i])
    
    z = float(hdu[0].header['RSHIFT'])
    wave = hdu[1].data['LAMBDA']/(1+z)
#    if wave[0]< 3395.:
#        continue
    flux = hdu[1].data['FLUX']
    err = np.sqrt(1/hdu[1].data['IVAR'])
    
    if 0.12<=z<=1.42:
        snr,peak = fit_line(fitslist[i][:-5],wave,flux,err,z)
        SNR.append(snr)
        Peak.append(peak)
        name.append(fitslist[i])
        
        
    else:
        continue
    print(i)
    
    
    
H = Table([name,SNR,Peak],names=['name','SNR','Peak'])
H.write('/Users/k/Documents/MMT_spec/MMT_LockmanHole/NeV/snr_table.ipac',format='ipac')
        
    
    
    
