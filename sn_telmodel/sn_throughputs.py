#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:50:53 2024

@author: philippe.gris@clermont.in2p3.fr
"""
from functools import wraps
from scipy.constants import *
from rubin_sim.phot_utils import photometric_parameters
from rubin_sim.phot_utils import Bandpass, Sed
import numpy as np
from sn_telmodel.sn_telescope import Telescope
from sn_telmodel.sn_atmosphere import Atmos_Transmission

# decorator to access parameters of the class


def get_val_decor(func):
    @wraps(func)
    def func_deco(theclass, what, xlist):
        for x in xlist:
            if x not in theclass.data[what].keys():
                func(theclass, what, x)
    return func_deco


def get_val_decorb(func):
    @wraps(func)
    def func_decob(theclass, what, xlist, y):
        for x in xlist:
            if x not in theclass.data[what].keys():
                func(theclass, what, x, y)
    return func_decob


class Throughputs(Telescope, Atmos_Transmission):
    def __init__(self, tel_dir='throughputs_1.9/baseline',
                 tel_optical_files=['detector.dat', 'lens1.dat',
                                    'lens2.dat', 'lens3.dat',
                                    'm1.dat', 'm2.dat', 'm3.dat'],
                 tel_filter_files=['filter_u.dat', 'filter_g.dat',
                                   'filter_r.dat', 'filter_i.dat',
                                   'filter_z.dat', 'filter_y.dat'],
                 tel_wave_min=300.,
                 tel_wave_max=1150.,
                 filter_colors=dict(zip(['u', 'g', 'r', 'i', 'z', 'y'],
                                        ['b', 'c', 'g', 'y', 'r', 'm'])),
                 site_name='LSST', pressure=743.,
                 atmos_dir='throughputs_1.9/atmos',
                 atmos_type='obsatmo',
                 darksky_file='throughputs_1.9/baseline/darksky.dat',
                 gain=2.5):
        """
        Throughputs class - inheritance from Telescope and Atmos_Transmission

        Parameters
        ----------
        tel_dir : str, optional
            input files directory. The default is 'throughputs_1.9/baseline'.
        tel_optical_files : list(str), optional
            List of optical components. The default is
                                        ['detector.dat', 'lens1.dat',
                                         'lens2.dat', 'lens3.dat',
                                         'm1.dat', 'm2.dat', 'm3.dat'].
        tel_filter_files : list(str), optional
            List of filter files. The default is
                                        ['filter_u.dat', 'filter_g.dat',
                                         'filter_r.dat', 'filter_i.dat',
                                        'filter_z.dat', 'filter_y.dat'].
        tel_wave_min : float, optional
            Min wave length. The default is 300..
        tel_wave_max : float, optional
            Max wavelength. The default is 1150..
        filter_colors : list(str), optional
            List of filter colors (for display). The default is
                dict(zip(['u', 'g', 'r', 'i', 'z', 'y'],                                     
                         ['b', 'c', 'g', 'y', 'r', 'm'])).
        site_name : str, optional
             Site name. The default is 'LSST'.
        pressure : float, optional
             Pressure. The default is 743..
        atmos_dir : str, optional
             Location dir of atmosphere files.
             The default is 'throughput_v1.9/atmos'.
        atmos_type : str, optional
             Transmission estimation method . The default is 'obsatmo'. 
        darksky_file: str, optional
            dark sky file. The default is 'throughputs_1.9/baseline/darksky.dat'
        gain: float, optional
            electronic gain. The default is 2.5.

        Returns
        -------
        None.

        """
        Telescope.__init__(self, tel_dir, tel_optical_files,
                           tel_filter_files, tel_wave_min, tel_wave_max,
                           filter_colors)
        Atmos_Transmission.__init__(self, site_name, pressure, atmos_dir,
                                    atmos_type)

        # load darksky
        self.load_darksky(darksky_file)

        # load_atmosphere
        self.load_atmosphere()

        # get throughputs
        self.throughputs = self.get_throughputs(self.atmosphere)

        # throughgputs data
        self.data = {}
        self.params = ['mag_sky', 'm5', 'Tb',
                       'Sigmab', 'zp', 'counts_zp', 'adu_zp',
                       'Skyb', 'flux_sky']
        self.reset_data()

        # electronic gain
        self.gain = gain

        # seeing
        self.data['FWHMeff'] = dict(
            zip('ugrizy', [0.92, 0.87, 0.83, 0.80, 0.78, 0.76]))

        # mean wavelength filters
        self.mean_wavelength = {}

    def reset_data(self):
        """
        Method to reset throughputs (zp, etc) data

        Returns
        -------
        None.

        """

        for par in self.params:
            self.data[par] = {}

    def new_atmosphere(self, site_name='LSST', airmass=1.2, aerosol=0.0,
                       pwv=4.0, oz=300, beta=1.4, pressure=743.):
        """
        Method to load atmospheric transmission

        Parameters
        ----------
        site_name : str, optional
            Site name. The default is 'LSST'.
        airmass : float, optional
            airmass. The default is 1.2.
        aerosol : float, optional
            aerosol. The default is 0.0.
        pwv : float, optional
            precipitable water vapor. The default is 4.0.
        oz : float, optional
            ozone. The default is 300.
        beta : float, optional
            Angstrom exponent. The default is 1.4.
        pressure : float, optional
            pressure. The default is 743..

        Returns
        -------
        None.

        """

        # load new atmosphere (update self.atmosphere)
        self.load_atmosphere(site_name, airmass, aerosol,
                             pwv, oz, beta, pressure)

        # get new throughputs
        self.throughputs = self.get_throughputs(self.atmosphere)

    def get_throughputs(self, bandpass):
        """
        Method to get the resulting throughputs=system*bandpass

        Parameters
        ----------
        bandpass : Bandpass
            Wavelength and sb.

        Returns
        -------
        through : dict
            Resulting troughput (system*filter).

        """

        through = {}
        for f in self.filter_list:
            wavelen, sb = self.tel_trans[f].multiply_throughputs(
                bandpass.wavelen, bandpass.sb)
            through[f] = Bandpass(wavelen=wavelen, sb=sb)

        return through

    def plot_throughputs(self, plt, fig=None, ax=None):
        """
        To plot the throughputs

        Parameters
        ----------
        plt : matplotlib.pyplot
            plot lib.
        fig : matplotlib figure, optional
            Figure for the plot. The default is None.
        ax : matplotlib axis, optional
            Axis for the plot. The default is None.

        Returns
        -------
        None.

        """

        # colors=['b','g','r','m','c',[0.8,0,0]]
        # style = [',', ',', ',', ',']

        if fig is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        for i, band in enumerate(self.filter_list):

            ax.plot(self.tel_trans[band].wavelen,
                    self.tel_trans[band].sb,
                    linestyle='--', color=self.filter_colors[band],
                    label='%s - tel' % (band))

            ax.plot(self.throughputs[band].wavelen,
                    self.throughputs[band].sb,
                    linestyle='-',
                    color=self.filter_colors[band],
                    label='%s - tel+atmos' % (band))

        # ax.plot(self.atmos.wavelen, self.atmos.sb, color='k',
        #        label='X =%.1f atmos' % (self.airmass), linestyle='-')

        ax.plot(self.atmosphere.wavelen, self.atmosphere.sb,
                color='k',
                label='atmos (airmass={})'.format(self.airmass),
                linestyle='--')
        # plt.legend(loc=(0.85, 0.1), fontsize='smaller',
        # fancybox=True, numpoints=1)

        ax.legend(loc=(0.82, 0.1), fancybox=True, numpoints=1)

        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Sb (0-1)')
        ax.set_title('System throughput')
        ax.grid(visible=True)

    def load_darksky(self, darksky_file):
        """
        Method to load the dark sky file

        Parameters
        ----------
        darksky_file : str
            dark sky file.

        Returns
        -------
        None.

        """

        self.darksky = Sed()
        self.darksky.read_sed_flambda(darksky_file)

    def plot_darksky(self, plt, fig=None, ax=None):
        """
        Method to plot the dark sky sed

        Parameters
        ----------
        plt : matplotlib pyplot
            plot lib.
        fig : matplotlib figure, optional
            plot figure. The default is None.
        ax : matplotlib axis, optional
            plot axis. The default is None.

        Returns
        -------
        None.

        """

        if fig is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        ax.plot(self.darksky.wavelen,
                self.darksky.flambda, 'k:', linestyle='-')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('flambda (ergs/cm$^2$/s/nm)')
        fig.suptitle('Dark Sky SED')
        ax.grid(visible=True)

    @get_val_decorb
    def get(self, what, band, exptime):
        """
        Decorator to access quantities

        Parameters
        ---------------
        what: str
          parameter to estimate
        band: str
          filter

        """
        filter_trans = self.throughputs[band]
        # wavelen_min, wavelen_max, wavelen_step = \
        #    filter_trans.get_wavelen_limits(None, None, None)
        wavelen_min = np.min(filter_trans.wavelen)
        wavelen_max = np.max(filter_trans.wavelen)
        wavelen_step = filter_trans.wavelen[1] - filter_trans.wavelen[0]

        # bpass = Bandpass(wavelen=filter_trans.wavelen, sb=filter_trans.sb)
        """
        flatSedb = Sed()
        flatSedb.set_flat_sed(wavelen_min, wavelen_max, wavelen_step)
        flux0b = 10.**(-0.4*(self.mag_sky(band)))
        flux0b = flatSedb.calc_flux_norm(self.mag_sky(band), filter_trans)
        flatSedb.multiply_flux_norm(flux0b)
        photParams = photometric_parameters.PhotometricParameters(
            bandpass=band)
        trans = filter_trans
        adu_int = flatSedb.calc_adu(bandpass=trans, phot_params=photParams)

        vv = adu_int
        vv *= photParams.platescale**2*photParams.gain
        vv /= (photParams.exptime)
        """
        photParams = photometric_parameters.PhotometricParameters(
            gain=self.gain, bandpass=band)
        photParams._exptime = exptime
        photParams._nexp = 1
        exptime = photParams.exptime
        nexp = photParams.nexp
        vv = self.mag_to_flux_e_sec(self.mag_sky(band), band, exptime, nexp)
        vv = vv[1]*photParams.platescale**2
        vvb = 10**(-0.4*(self.mag_sky(band)-self.zp(band)))
        vvb *= photParams.platescale**2

        print(vv, vvb)
        self.data['flux_sky'][band] = vv

        trans = self.throughputs[band]

        from rubin_sim.phot_utils import signaltonoise
        nexp = exptime/30
        photParams._nexp = nexp
        photParams._exptime = exptime/nexp

        flatSedb = Sed()
        flatSedb.set_flat_sed(wavelen_min, wavelen_max, wavelen_step)
        flux0b = flatSedb.calc_flux_norm(self.mag_sky(band), filter_trans)
        flatSedb.multiply_flux_norm(flux0b)

        self.data['m5'][band] = signaltonoise.calc_m5(
            flatSedb, trans, filter_trans,
            phot_params=photParams,
            fwhm_eff=self.FWHMeff(band))

    @get_val_decor
    def get_inputs(self, what, band):
        """
        decorator to access Tb, Sigmab, mag_sky

        Parameters
        ---------------
        what: str
          parameter to estimate
        band: str
          filter

        """
        myup = self.Calc_Integ_Sed(self.darksky, self.throughputs[band])
        # bpass = self.atmosphere[band]
        # if self.aerosol_b:
        bpass = self.throughputs[band]
        self.data['Tb'][band] = self.Calc_Integ(bpass)
        self.data['Sigmab'][band] = self.Calc_Integ(self.throughputs[band])
        tt = np.log10(myup/(3631.*self.Sigmab(band)))
        self.data['mag_sky'][band] = -2.5 * tt

    @get_val_decor
    def get_zp(self, what, band):
        """
        decorator get zero points
        formula used here are extracted from LSE-40

        Parameters
        ---------------
        band: str
          filter

        """
        photParams = photometric_parameters.PhotometricParameters(gain=self.gain,
                                                                  bandpass=band)
        photParams._exptime = 30
        photParams._nexp = 1
        Diameter = 2.*np.sqrt(photParams.effarea*1.e-4 /
                              np.pi)  # diameter in meter
        Cte = 3631.*np.pi*Diameter**2*photParams.exptime/4/h/1.e36

        self.data['Skyb'][band] = Cte*np.power(Diameter/6.5, 2.)\
            * np.power(photParams.exptime/30., 1.)\
            * np.power(photParams.platescale, 2.)\
            * 10.**0.4*(25.-self.mag_sky(band))\
            * self.Sigmab(band)

        Zb = 181.8*np.power(Diameter/6.5, 2.)*self.Tb(band)
        mbZ = 25.+2.5*np.log10(Zb)
        self.data['zp'][band] = mbZ

        """
        filtre_trans = self.atmosphere[band]
        if self.atmos:
            filtre_trans = self.atmosphere[band]
        if self.aerosol_b:
            filtre_trans = self.aerosol[band]
        """
        filtre_trans = self.throughputs[band]
        # filtre_trans = self.aerosol[band]
        """
        wavelen_min, wavelen_max, wavelen_step = \
            filtre_trans.get_wavelen_limits(None, None, None)

        #filtre_trans = self.lsst_system[band]
        wavelen_min, wavelen_max, wavelen_step = filtre_trans.get_wavelen_limits(
            None, None, None)
        """
        wavelen_min = np.min(filtre_trans.wavelen)
        wavelen_max = np.max(filtre_trans.wavelen)
        wavelen_step = filtre_trans.wavelen[1] - filtre_trans.wavelen[0]

        # diff = np.diff(filtre_trans.wavelen)
        bpass = Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
        flatSed = Sed()
        flatSed.set_flat_sed(wavelen_min,
                             wavelen_max, wavelen_step)

        flux0 = np.power(10., -0.4*mbZ)
        flatSed.multiply_flux_norm(flux0)
        photParams = photometric_parameters.PhotometricParameters(
            gain=self.gain, bandpass=band)

        # number of counts for exptime
        counts = flatSed.calc_adu(bpass, phot_params=photParams)

        self.data['counts_zp'][band] = counts*photParams.gain / \
            (photParams.exptime*photParams.nexp)

        self.data['adu_zp'][band] = counts / \
            (photParams.exptime*photParams.nexp)

    def return_value(self, what, band):
        """
        accessor

        Parameters
        ---------------
        what: str
          parameter to estimate
        band: str
          filter

        """
        if len(band) > 1:
            return self.data[what]
        else:
            return self.data[what][band]

    def m5(self, filtre, exptime):
        """m5 accessor
        """
        self.get('m5', filtre, exptime)
        return self.return_value('m5', filtre)

    def Tb(self, filtre):
        """Tb accessor
        """
        self.get_inputs('Tb', filtre)
        return self.return_value('Tb', filtre)

    def mag_sky(self, filtre):
        """mag_sky accessor
        """
        self.get_inputs('mag_sky', filtre)
        return self.return_value('mag_sky', filtre)

    def flux_sky(self, filtre, exptime):
        """flux_sky accessor
        """

        self.get('flux_sky', filtre, exptime)
        return self.return_value('flux_sky', filtre)

    def Sigmab(self, filtre):
        """
        Sigmab accessor

        Parameters
        ----------------
        band: str
          filter

        """
        self.get_inputs('Sigmab', filtre)
        return self.return_value('Sigmab', filtre)

    def zp(self, filtre):
        """
        zp accessor

        Parameters
        ----------------
        band: str
          filter

        """
        self.get_zp('zp', filtre)
        return self.return_value('zp', filtre)

    def counts_zp(self, filtre):
        """
        counts_zp accessor

        Parameters
        ----------
        filtre : str
            filter to consider.

        Returns
        -------
        None.

        """
        self.get_zp('zp', filtre)
        return self.return_value('counts_zp', filtre)

    def adu_zp(self, filtre):
        """
        counts_zp accessor

        Parameters
        ----------
        filtre : str
            filter to consider.

        Returns
        -------
        None.

        """
        self.get_zp('zp', filtre)
        return self.return_value('adu_zp', filtre)

    def FWHMeff(self, filtre):
        """
        FWHMeff accessor

        Parameters
        ----------------
        band: str
          filter
        """
        return self.return_value('FWHMeff', filtre)

    def Calc_Integ(self, bandpass):
        """
        integration over bandpass

        Parameters
        --------------
        bandpass : float

        Returns
        ---------
        integration

        """
        resu = 0.
        dlam = 0
        for i, wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam = bandpass.wavelen[i+1]-wave
                resu += dlam*bandpass.sb[i]/wave
            # resu+=dlam*bandpass.sb[i]

        return resu

    def Calc_Integ_Sed(self, sed, bandpass, wavelen=None, fnu=None):
        """
        SED integration

        Parameters
        --------------
        sed : float
          sed to integrate
        bandpass : float
          bandpass
        wavelength : float, opt
          wavelength values
           Default : None
        fnu : float, opt
           fnu values
           Default : None

        Returns
        ----------
        integrated sed over the bandpass

        """
        use_self = sed._check_use_self(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if sed.fnu is None:
                # If fnu not present, calculate. (does not regrid).
                sed.flambda_tofnu()
            wavelen = sed.wavelen
            fnu = sed.fnu
        # Make sure wavelen/fnu are on the same wavelength grid as bandpass.
        wavelen, fnu = sed.resample_sed(
            wavelen, fnu, wavelen_match=bandpass.wavelen)

        # Calculate the number of photons.
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        dlambda = wavelen[1] - wavelen[0]
        return nphoton * dlambda

    def flux_to_mag(self, flux, band, zp=None):
        """
        Flux to magnitude conversion

        Parameters
        --------------
        flux : float
          input fluxes
        band : str
           input band
        zp : float, opt
           zeropoints
           Default : None

        Returns
        ---------
        magnitudes

        """
        if zp is None:
            zp = self.zero_points(band)
        # print 'zp',zp,band
        m = -2.5 * np.log10(flux) + zp
        return m

    def mag_to_flux(self, mag, band, zp=None):
        """
        Magnitude to flux conversion

        Parameters
        --------------
        mag : float
          input mags
        band : str
           input band
        zp : float, opt
           zeropoints
           Default : None

        Returns
        ---------
        fluxes

        """
        if zp is None:
            zp = self.zero_points(band)
        return np.power(10., -0.4 * (mag-zp.item()))

    def zero_points(self, band):
        """
        Zero points estimation

        Parameters
        --------------
        band : list(str)
          list of bands

        Returns
        ---------
        array of zp

        """
        return np.asarray([self.zp(b) for b in band])

    def mag_to_flux_e_sec(self, mag, band, exptime, nexp):
        """
        Mag to flux (in photoelec/sec) conversion

        Parameters
        --------------
        mag : float
          input magnitudes
        band : str
          input bands
        exptime : float
          input exposure times
        nexp: int
          number of exposures

        Returns
        ----------
        counts : float
           number of ADU counts
        e_per_sec : float
           flux in photoelectron per sec.

        """
        filter_trans = self.throughputs[band]
        if not hasattr(mag, '__iter__'):

            # wavelen_min, wavelen_max, wavelen_step = filter_trans.get_wavelen_limits(
            #    None, None, None)
            wavelen_min = np.min(filter_trans.wavelen)
            wavelen_max = np.max(filter_trans.wavelen)
            wavelen_step = filter_trans.wavelen[1] - filter_trans.wavelen[0]

            sed = Sed()
            sed.set_flat_sed()

            flux0 = sed.calc_flux_norm(mag, filter_trans)
            sed.multiply_flux_norm(flux0)

            photParams = photometric_parameters.PhotometricParameters(gain=self.gain,
                                                                      exptime=exptime, nexp=nexp)

            counts = sed.calc_adu(
                bandpass=filter_trans, phot_params=photParams)
            e_per_sec = counts

            # counts per sec
            e_per_sec /= exptime
            counts /= exptime
            # conversion to pe
            e_per_sec *= photParams.gain
            return counts, e_per_sec
        else:
            r = []
            for m, b, expt, nexpos in zip(mag, band, exptime, nexp):
                counts, flux_e = self.mag_to_flux_e_sec(m, b, expt, nexpos)
                r.append((counts, flux_e))
            return np.asarray(r)

    def gamma(self, mag, band, exptime, nexp):
        """
        gamma parameter estimation

        cf eq(5) of the paper LSST :
            from science drivers to reference design
            and anticipated data products

        with sigma_rand = 0.2 and m=m5

        Parameters
        --------------
        mag : float
          magnitudes
        band : str
          band
        exptime : float
          exposure time

        Returns
        ----------
        gamma, mag_to_flux (float)

        """

        if not hasattr(mag, '__iter__'):
            photParams = photometric_parameters.PhotometricParameters(gain=self.gain,
                                                                      nexp=nexp, exptime=exptime)
            counts, e_per_sec = self.mag_to_flux_e_sec(
                mag, band, exptime, nexp)
            gamma = 0.04-1./(photParams.gain*counts)
            return gamma, e_per_sec
        else:
            r = []
            for m, b, e, nexpo in zip(mag, band, exptime, nexp):
                gamma, flux_e = self.gamma(m, b, e, nexpo)
                r.append((gamma, flux_e))
            return np.asarray(r)

    def etc(self, exptime=30., plateScale=0.2):
        """
        Method to print the throughputs parameters

        Parameters
        ----------
        exptime : float, optional
            exposure time. The default is 30..
        plateScale : float, optional
            plate scale ("2). The default is 0.2.

        Returns
        -------
        None.

        """

        import pandas as pd
        # exptime = 30
        # plateScale = 0.2  # pixel size ''
        bands = self.filter_list
        df = pd.DataFrame(list(bands), columns=['band'])
        zp = dict(zip(bands, [self.zp(b) for b in bands]))
        mag_sky = dict(zip(bands, [self.mag_sky(b) for b in bands]))
        flux_sky = dict(zip(bands, [self.flux_sky(b, exptime) for b in bands]))
        m5 = dict(zip(bands, [self.m5(b, exptime) for b in bands]))

        df['zp'] = [self.zp(b) for b in bands]
        df['flux_zp'] = [self.counts_zp(b) for b in bands]
        df['ADU_zp'] = [self.adu_zp(b) for b in bands]
        df['msky'] = [self.mag_sky(b) for b in bands]
        df['flux_sky'] = [self.flux_sky(b, exptime) for b in bands]
        df['flux_sky_from_mag'] = 10**(-0.4 *
                                       (df['msky']-df['zp']))*plateScale**2
        # df['flux_sky_mag'] = -2.5*np.log10(df['flux_sky'])+df['zp']
        df['FWHMeff'] = [self.FWHMeff(b) for b in bands]
        df['m5'] = [self.m5(b, exptime) for b in bands]

        df = df.rename(columns={"zp": "zp (AB)",
                                "flux_zp": "flux_zp (pe/s/pix)",
                                "flux_sky": "flux_sky (pe/s/pix)",
                                "flux_sky_from_mag": "flux_sky_from_mag (pe/s/pix)",
                                "FWHMeff": "FWHMEff ('')",
                                "m5": "m5 (exptime: {} s)".format(exptime),
                                "msky": "msky (/\"2)"})
        df = df.round(2)
        pd.set_option('display.colheader_justify', 'center')
        print(df.to_string(index=False))

    def mean_wave(self):
        """ Estimate mean wave
        """
        for band in self.filter_list:
            self.mean_wavelength[band] = np.sum(
                self.throughputs[band].wavelen*self.throughputs[band].sb)\
                / np.sum(self.throughputs[band].sb)


def load_throughputs_from_config(config):
    """
    Function to load telescope model

    Parameters
    ----------
    config : dict
        configuration parameters.

    Returns
    -------
    tel : Telescope class (sn_telmodel.sn_telescope)
        Telescope model.

    """

    name = config['name']
    tel_dir = config['telescope']['dir']
    tel_tag = config['telescope']['tag']
    through_dir = config['throughputDir']
    atmos_dir = config['atmosDir']
    airmass = config['airmass']
    aerosol = config['aerosol']
    pwv = config['pwv']
    oz = config['oz']

    # point_to_tag(tel_dir, tel_tag)

    tel_dir = '{}_{}'.format(tel_dir, tel_tag)
    through_dir = '{}/{}'.format(tel_dir, through_dir)
    atmos_dir = '{}/{}'.format(tel_dir, atmos_dir)

    airmass = float(airmass)
    aerosol = float(aerosol)
    pwv = float(pwv)
    oz = float(oz)

    throughputs = Throughputs(tel_dir=through_dir,
                              site_name=name,
                              atmos_dir=atmos_dir,
                              atmos_type='obsatmo')

    throughputs.new_atmosphere(site_name=name,
                               airmass=airmass,
                               aerosol=aerosol,
                               pwv=pwv, oz=oz)

    """
    tel = get_telescope(name=name, tel_dir=tel_dir,
                        through_dir=through_dir,
                        atmos_dir=atmos_dir, airmass=airmass,
                        aerosol=aerosol, pwv=pwv, oz=oz, tag=tel_tag)
    """
    return throughputs


def get_telescope(name='LSST',
                  tel_dir='throughputs',
                  through_dir='baseline',
                  atmos_dir='atmos',
                  tag='1.9', airmass=1.2, gain=2.5,
                  pwv=4.0, oz=400,
                  aerosol=0.0, beta=1.4, pressure=743.,
                  load_components=False):
    """
    Function to grab telescope version

    Parameters
    ----------
    name : str, optional
       Telescope name. The default is 'LSST'.
    tel_dir : str, optional
       Main tel directory. The default is 'throughputs'.
    through_dir : str, optional
        Throughput directory. The default is 'throughputs/baseline'.
    atmos_dir : str, optional
        Atmosphere directory. The default is 'throughputs/atmos'.
    tag : str, optional
        Tag version for throughputs. The default is '1.9'.
    airmass : float, optional
        airmass value for throughputs. The default is 1.2.
    gain: float, optional.
         electronic gain. The default is 2.5
    aerosol : bool, optional
        add aerosol effect. The default is True.
    load_components : bool, optional
        To load all the components (one by one). The default is False.
    Returns
    -------
    tela : TYPE
        DESCRIPTION.

    """

    # print('Telescope instance', tel_dir)

    """
    tel = Telescope(tel_dir=tel_dir,
                    airmass=airmass, through_dir=through_dir,
                    atmos_dir=atmos_dir, aerosol=aerosol, pwv=pwv, oz=oz,
                    beta=beta, pressure=pressure,
                    load_components=load_components, tag=tag, gain=gain)
    """
    """
    tel = Throughputs(tel_dir=tel_dir,
                      airmass=airmass, through_dir=through_dir,
                      atmos_dir=atmos_dir, aerosol=aerosol, pwv=pwv, oz=oz,
                      beta=beta, pressure=pressure,
                      load_components=load_components, tag=tag, gain=gain)
    """

    throughputs = Throughputs(tel_dir=through_dir,
                              site_name=name,
                              atmos_dir=atmos_dir,
                              atmos_type='obsatmo')

    throughputs.new_atmosphere(site_name=name,
                               airmass=airmass,
                               aerosol=aerosol,
                               pwv=pwv, oz=oz)

    return throughputs


def point_to_tag(tel_dir, tag):
    """
    Function to point to a given tel tag version

    Parameters
    ----------
    tel_dir : str
        Main telescope dir.
    tag : str
        Tag throughputs version.

    Returns
    -------
    None.

    """

    import os
    path = os.getcwd()
    throughputs_dir = '{}_{}'.format(tel_dir, tag)
    if not os.path.isdir(throughputs_dir):
        cmd = 'git clone https://github.com/lsst/{} {}_{}'.format(
            tel_dir, tel_dir, tag)
        os.system(cmd)

        os.chdir(throughputs_dir)
        cmd = 'git checkout tags/{}'.format(tag)
        os.system(cmd)
        os.chdir(path)
