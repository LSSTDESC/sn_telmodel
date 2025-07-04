#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:24:45 2024

@author: philippe.gris@clermont.in2p3.fr
"""

import numpy as np
import pandas as pd


def get_trans(am, pwv, oz, tau=0., beta=1.4,
              colname=['Wavelength(nm)', 'Throughput(0-1)'], emul=None):
    """
    Function to estimate transmissions

    Parameters
    ----------
    am : float
        airmass value.
    pwv : float
        precipitable water vapor in mm.
    oz : float
        Ozone depth in DU (Dobson Unit).
    tau : float, optional
        vertical aerosol depth of each component at lambda0 vavelength.
        The default is 0..
    beta : float, optional
       the angstrom exponent. Must be positive in the range 0., 3.
       The default is 1.4.
    colname : list(str), optional
        list of output columns. 
        The default is ['Wavelength(nm)', 'Throughput(0-1)'].

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    # emulate obsAtmo
    if emul is None:
        from getObsAtmo.getObsAtmo import ObsAtmo
        emul = ObsAtmo('LSST', 743.0)
    wl = list(np.arange(300., 1100., 0.1))
    transm = emul.GetAllTransparencies(wl, am, pwv, oz, tau=tau, beta=beta)
    df = pd.DataFrame(wl, columns=[colname[0]])
    df[colname[1]] = transm

    """
    # extend down to 200 nm
    wl = list(np.arange(250., 300., 0.5))
    dfb = pd.DataFrame(wl, columns=[colname[0]])
    dfb[colname[1]] = 0.0
    """
    # extend up  to 1150 nm
    wl = list(np.arange(1100, 1150.1, 0.1))
    dfc = pd.DataFrame(wl, columns=[colname[0]])
    idx = np.abs(df[colname[0]]-1100.) < 0.01
    dfc[colname[1]] = df[colname[1]].iloc[-1]

    # df = pd.concat((df, dfb))
    df = pd.concat((df, dfc))

    df = df.sort_values(by=colname[0])
    return df


class Zeropoint_airmass:
    def __init__(self, throughputs, pwv=4.0, ozone=400., aerosol=0.0):
        """
        class to estimate zp vs airmass and fit (linear) the results

        Parameters
        ----------
        throughputs: Throughputs class
          instance of a throughput class
        pwv : float, optional
            pwv value. The default is 4.0.
        ozone :float, optional
            ozone value. The default is 400..
        aerosol : float, optional
            aerosol value. The default is 0.0.
        Returns
        -------
        None.

        """

        self.throughputs = throughputs
        self.pwv = pwv
        self.ozone = ozone
        self.aerosol = aerosol

    def get_data(self):
        """
        Method to estimate zp vs airmass

        Returns
        -------
        res : numpy array
            columns: filter, zp, airmass, mean_wave.

        """

        r = []
        # point_to_tag(self.tel_dir, self.tag)
        tel = self.throughputs
        for airmass in np.arange(1., 2.51, 0.1):
            """
            tel = get_telescope(tel_dir=tel_dir,
                                through_dir=through_dir,
                                atmos_dir=atmos_dir,
                                tag=self.tag, load_components=True,
                                airmass=airmass,
                                aerosol=self.aerosol, pwv=self.pwv, oz=self.oz)
            """
            tel.new_atmosphere(site_name=tel.site_name,
                               airmass=airmass,
                               aerosol=self.aerosol,
                               pwv=self.pwv, ozone=self.ozone)
            tel.mean_wave()
            for b in 'ugrizy':
                # b = 'g'
                # print(airmass, b, tel.zp(b))
                mean_wave = tel.mean_wavelength[b]
                rb = [airmass]
                rb.append(b)
                rb.append(tel.zp(b))
                rb.append(tel.counts_zp(b))
                rb.append(mean_wave)
                r.append(rb)

        res = np.rec.fromrecords(
            r, names=['airmass', 'band', 'zp', 'zp_e_sec', 'mean_wavelength'])

        return res

    def fitfunc(self, x, a, b):
        """
        Function used for fitting

        Parameters
        ----------
        x : array(float)
            x-axis var.
        a : float
            slope.
        b : float
            intercept.

        Returns
        -------
        array
            list of values.

        """

        return a*x+b

    def fit(self, res, xvar='airmass', yvar='zp'):
        """
        Function to fit yvar vs xvar for all bands.

        Parameters
        ----------
        res : array
            data to fit.
        xvar : str, optional
            x-axis var. The default is 'airmass'.
        yvar : str, optional
            y-axis var. The default is 'zp'.

        Returns
        -------
        res : array
            slop and intercep from the fit per band.
            added mean_wavelength.

        """

        from scipy.optimize import curve_fit
        r = []
        for b in 'ugrizy':
            idx = res['band'] == b
            sel = res[idx]
            xdata = sel[xvar]
            ydata = sel[yvar]
            popt, pcov = curve_fit(self.fitfunc, xdata, ydata)
            mean_wave = np.mean(sel['mean_wavelength'])
            r.append((b, popt[0], popt[1], mean_wave))

        res = np.rec.fromrecords(
            r, names=['band', 'slope', 'intercept', 'mean_wavelength'])

        return res

    def get_fit_params(self):

        # get data
        data = self.get_data()

        # fit these data
        fitdata = self.fit(data)

        return fitdata
