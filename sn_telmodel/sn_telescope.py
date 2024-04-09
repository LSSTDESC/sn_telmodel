from functools import wraps
from scipy.constants import *
import numpy as np
from rubin_sim.phot_utils import photometric_parameters
from rubin_sim.phot_utils import Bandpass, Sed
from sn_telmodel.sn_throughputs import Throughputs
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# import matplotlib.pyplot as plt
# import math

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


class Zeropoint_airmass:
    def __init__(self, tel_dir='throughputs',
                 through_dir='baseline',
                 atmos_dir='atmos',
                 tag='1.9', aerosol=True):
        """
        class to estimate zp vs airmass and fit (linear) the results

        Parameters
        ----------
        tel_dir : str, optional
            telescope main dir. The default is 'throughputs'.
        through_dir : str, optional
            throughputs dir. The default is 'throughputs/baseline'.
        atmos_dir : str, optional
            dir for atmos files. The default is 'throughputs/atmos'.
        tag : str, optional
            throughputs tag. The default is '1.9'.
        aerosol : bool, optional
            to include aerosol. The default is True.

        Returns
        -------
        None.

        """

        self.tel_dir = tel_dir
        self.through_dir = through_dir
        self.atmos_dir = atmos_dir
        self.tag = tag
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
        point_to_tag(self.tel_dir, self.tag)
        tel_dir = '{}_{}'.format(self.tel_dir, self.tag)
        through_dir = '{}/{}'.format(tel_dir, self.through_dir)
        atmos_dir = '{}/{}'.format(tel_dir, self.atmos_dir)
        for airmass in np.arange(1., 2.51, 0.1):
            tel = get_telescope(tel_dir=tel_dir,
                                through_dir=through_dir,
                                atmos_dir=atmos_dir,
                                tag=self.tag, airmass=airmass,
                                aerosol=self.aerosol)
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
            r, names=['airmass', 'band', 'zp', 'zp_adu_sec', 'mean_wavelength'])

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


def load_telescope_from_config(config):
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

    point_to_tag(tel_dir, tel_tag)

    tel_dir = '{}_{}'.format(tel_dir, tel_tag)
    through_dir = '{}/{}'.format(tel_dir, through_dir)
    atmos_dir = '{}/{}'.format(tel_dir, atmos_dir)

    tel = get_telescope(name=name, tel_dir=tel_dir,
                        through_dir=through_dir,
                        atmos_dir=atmos_dir, airmass=airmass,
                        aerosol=aerosol, tag=tel_tag)

    return tel


def get_telescope(name='LSST',
                  tel_dir='throughputs',
                  through_dir='baseline',
                  atmos_dir='atmos',
                  tag='1.9', airmass=1.2, gain=2.5,
                  aerosol=True, load_components=False):
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

    print('Telescope instance', tel_dir)

    tel = Telescope(name=name, tel_dir=tel_dir,
                    airmass=airmass, through_dir=through_dir,
                    atmos_dir=atmos_dir, aerosol=aerosol,
                    load_components=load_components, tag=tag, gain=gain)

    return tel


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


class Telescope(Throughputs):
    """ Telescope class
    inherits from Throughputs
    estimate quantities defined in LSE-40

    The following quantities are accessible:

    mag_sky: sky magnitude

    m5: 5-sigma depth

    Sigmab: see eq. (36) of LSE-40

    zp: see eq. (43) of LSE-40

    counts_zp:

    Skyb: see eq. (40) of LSE-40

    flux_sky:


    Parameters
    -------------
    through_dir : str, opt
       throughput directory
       Default : LSST_THROUGHPUTS_BASELINE
    atmos_dir : str, opt
       directory of atmos files
       Default : THROUGHPUTS_DIR
    telescope_files : list(str),opt
       list of of throughput files
       Default : ['detector.dat', 'lens1.dat','lens2.dat',
           'lens3.dat','m1.dat', 'm2.dat', 'm3.dat']
    filterlist: list(str), opt
       list of filters to consider
       Default : 'ugrizy'
    wave_min : float, opt
        min wavelength for throughput
        Default : 300
    wave_max : float, opt
        max wavelength for throughput
        Default : 1150
    atmos : bool, opt
         to include atmosphere affects
         Default : True
    aerosol : bool, opt
         to include aerosol effects
         Default : True
    airmass : float, opt
         airmass value
         Default : 1.

    Returns
    ---------
    Accessible throughputs (per band, from Throughput class):
    lsst_system: system throughput (lens+mirrors+filters)
    lsst_atmos: lsst_system+atmosphere
    lsst_atmos_aerosol: lsst_system+atmosphere+aerosol


    """

    def __init__(self, name='unknown', airmass=1., aerosol='',
                 tel_dir='throughputs', tag='1.9', gain=2.5, **kwargs):
        super().__init__(**kwargs)
        """
        self.name = name

        Throughputs.__init__(self, **kwargs)
        """
        params = ['mag_sky', 'm5', 'FWHMeff', 'Tb',
                  'Sigmab', 'zp', 'counts_zp', 'Skyb', 'flux_sky']
        self.name = name
        self.tel_dir = tel_dir
        self.tag = tag
        self.data = {}
        self.gain = gain
        for par in params:
            self.data[par] = {}

        self.data['FWHMeff'] = dict(
            zip('ugrizy', [0.92, 0.87, 0.83, 0.80, 0.78, 0.76]))

        # self.data['FWHMeff'] = dict(
        #    zip('ugrizy', [0.77, 0.73, 0.70, 0.67, 0.65, 0.63]))

        # self.atmos = atmos

        self.load_atmosphere(airmass, aerosol)

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
        filter_trans = self.system[band]
        # filter_trans = self.lsst_atmos_aerosol[band]
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

        if self.atmos:
            trans = self.atmosphere[band]

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

    @ get_val_decor
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
        myup = self.Calc_Integ_Sed(self.darksky, self.system[band])
        bpass = self.atmosphere[band]
        if self.aerosol_b:
            bpass = self.aerosol[band]
        self.data['Tb'][band] = self.Calc_Integ(bpass)
        self.data['Sigmab'][band] = self.Calc_Integ(self.system[band])
        tt = np.log10(myup/(3631.*self.Sigmab(band)))
        self.data['mag_sky'][band] = -2.5 * tt

    @ get_val_decor
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

        filtre_trans = self.atmosphere[band]
        if self.atmos:
            filtre_trans = self.atmosphere[band]
        if self.aerosol_b:
            filtre_trans = self.aerosol[band]
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
        filter_trans = self.atmosphere[band]
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
