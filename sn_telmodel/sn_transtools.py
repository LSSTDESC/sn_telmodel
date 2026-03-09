#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 16:24:45 2024

@author: philippe.gris@clermont.in2p3.fr
"""

import numpy as np
import pandas as pd
from random import gauss


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
    def __init__(self, throughputs, pwv=4.0, ozone=400.,
                 aerosol=0.0, exptime=30., nexp=1):
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
        exptime: float, optional.
            exposure time [s]. the default is 30.
        nexp: float, optional
            number of exposure. The default is 1.
        Returns
        -------
        None.

        """

        self.throughputs = throughputs
        self.pwv = pwv
        self.ozone = ozone
        self.aerosol = aerosol
        self.exptime = exptime
        self.nexp = nexp

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
                rb.append(tel.zp(b, exptime=self.exptime, nexp=self.nexp))
                rb.append(tel.counts_zp(
                    b, exptime=self.exptime, nexp=self.nexp))
                rb.append(mean_wave)
                r.append(rb)
            tel.reset_data()

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


class Zeropoint_sigma_airmass:
    def __init__(self, config, exptime=30., nexp=1):
        """
        class to estimate zp,sigma zp vs airmass

        Parameters
        ----------
        config: dict
          config file for throughput params
        exptime: float, optional.
            exposure time [s]. the default is 30.
        nexp: float, optional
            number of exposure. The default is 1.
        Returns
        -------
        None.

        """
        from sn_telmodel.sn_throughputs import load_throughputs_from_config
        self.throughputs = load_throughputs_from_config(config)
        self.ntrial = config['ntrial']['zp']
        # airmass parameters
        self.airmass = config['airmass']
        self.sigma_airmass = config['sigma']['airmass']
        # airmass_round = config['round']['airmass']
        # pwv parameters
        self.pwv = config['pwv']
        self.sigma_pwv = config['sigma']['pwv']
        # pwv_round = config['round']['pwv']
        # ozone parameters
        self.ozone = config['ozone']
        self.sigma_ozone = config['sigma']['ozone']
        # ozone_round = config['round']['ozone']
        # aerosol parameters
        self.aerosol = config['aerosol']
        self.sigma_aerosol = config['sigma']['aerosol']

        self.exptime = exptime
        self.nexp = nexp

    def get_data(self):
        """
        Method to estimate zp, sigma_zp, mean_wave,sigma_mean_wave vs airmass

        Returns
        -------
        dict
            dict of interpolators.

        """

        from sn_tools.sn_utils import multiproc
        params = {}
        airmass = np.arange(1., 2.7, 0.1).tolist()

        import time
        time_ref = time.time()
        df = multiproc(airmass, params, self.get_param_loop, nproc=8)

        print('zeropoint+sigmas', time.time()-time_ref)

        """
        print(df)

        # cross check
        self.sigma_aerosol = 0
        self.sigma_pwv = 0
        self.sigma_ozone = 0.
        self.sigma_airmass = 0
        self.ntrial = 1

        dfi = self.get_data_indiv()

        print(dfi)
        print(test)
        """
        return self.interpIt(df)

    def interpIt(self, df):
        """
        Estimate ID interpolators

        Parameters
        ----------
        df : pandas df
            Data to make interp with.

        Returns
        -------
        dd : dict
            dict of interpolators.

        """

        from scipy.interpolate import interp1d

        dd = {}
        bands = df['band'].unique()
        for b in bands:
            idx = df['band'] == b
            sel = df[idx]

            for vv in ['zp', 'sigma_zp', 'mean_wave', 'sigma_mean_wave']:
                if vv not in dd.keys():
                    dd[vv] = {}
                dd[vv][b] = interp1d(sel['airmass'],
                                     sel[vv],
                                     bounds_error=False,
                                     fill_value=0.)
        return dd

    def get_param_loop(self, vals, params, j=0, output_q=None):
        """
         Method to estimate zp,mean_wave

         Parameters
         ----------
         vals : list(float)
             airmass values.
         params : dict
             parameters.
         j : int, optional
             tag for multiprocessing. The default is 0.
         output_q : multiprocessing queue, optional
             where to put the results. The default is None.

         Returns
         -------
         TYPE
             DESCRIPTION.

         """

        r = []
        # point_to_tag(self.tel_dir, self.tag)
        tel = self.throughputs
        df = pd.DataFrame()

        for airmass in vals:
            dfc = pd.DataFrame()

            for i in range(self.ntrial):
                dfa = self.get_params(tel, airmass)
                dfc = pd.concat((dfc, dfa))
                tel.reset_data
            dfb = dfc.groupby(['band']).apply(
                lambda x: self.stat(x),include_groups=False).reset_index()
            dfb['airmass'] = airmass
            df = pd.concat((df, dfb))

        if output_q is not None:
            return output_q.put({j: df})
        else:
            return df

    def get_data_indiv(self):
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
        df = pd.DataFrame()
        import time
        for airmass in np.arange(1., 2.81, 0.1):
            dfc = pd.DataFrame()
            time_ref = time.time()
            for i in range(self.ntrial):
                dfa = self.get_params(tel, airmass)
                dfc = pd.concat((dfc, dfa))
                tel.reset_data
            dfb = dfc.groupby(['band']).apply(
                lambda x: self.stat(x)).reset_index()
            dfb['airmass'] = airmass
            df = pd.concat((df, dfb))
            print('done', time.time()-time_ref)

        return df

    def get_params(self, tel, airmass):
        """
        Method to estimate zp, deltazp, mean_wave, delta_mean_wave for a given airmass value

        Parameters
        ----------
        tel : Throughput instance
            Throughputs to use.
        airmass : float
            airmass value.

        Returns
        -------
        res: pandas df
          output value

        """

        airmass += gauss(0, self.sigma_airmass)
        aerosol = self.aerosol+gauss(0, self.sigma_aerosol)
        pwv = self.pwv+gauss(0, self.sigma_pwv)
        ozone = self.ozone+gauss(0, self.sigma_ozone)

        if airmass < 1:
            return pd.DataFrame()
        tel.new_atmosphere(site_name=tel.site_name,
                           airmass=airmass,
                           aerosol=aerosol,
                           pwv=pwv, ozone=ozone)
        tel.mean_wave()
        r = []
        for b in 'ugrizy':
            # b = 'g'
            # print(airmass, b, tel.zp(b))
            mean_wave = tel.mean_wavelength[b]
            rb = [b]
            rb.append(tel.zp(b, exptime=self.exptime, nexp=self.nexp))
            rb.append(tel.counts_zp(
                b, exptime=self.exptime, nexp=self.nexp))
            rb.append(mean_wave)
            r.append(rb)
        tel.reset_data()

        res = pd.DataFrame(
            r, columns=['band', 'zp', 'zp_e_sec', 'mean_wave'])

        return res

    def stat(self, grp):
        """
        Method to estimate mean and rmses

        Parameters
        ----------
        grp : pandas df
            Data to process.

        Returns
        -------
        res : pandas df
            Mean dn rms.

        """

        dd = {}
        for vv in ['zp', 'mean_wave']:
            dd[vv] = [grp[vv].mean()]
            grp_std = 0.
            if len(grp) > 1:
                grp_std = grp[vv].std()
            dd['sigma_{}'.format(vv)] = [grp_std]

        res = pd.DataFrame.from_dict(dd)

        return res


def zp_from_config(config_instr):
    """
    Parameters
    ----------
    config_instr : dict
        config parameters.

    Returns
    -------
    interp
        zp vs airmass per band.

    """

    zp = Zeropoint_sigma_airmass(config_instr)

    return zp.get_data()
