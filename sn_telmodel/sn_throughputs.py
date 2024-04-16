import os
import matplotlib.pyplot as plt
from rubin_sim.phot_utils import Bandpass, Sed
# from rubin_sim.phot_utils import Sed
import numpy as np


class Throughputs(object):
    """ class to handle instrument throughput

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

    Returns
    ---------
    Accessible throughputs (per band):
    lsst_system: system throughput (lens+mirrors+filters)
    lsst_atmos: lsst_system+atmosphere
    lsst_atmos_aerosol: lsst_system+atmosphere+aerosol
    """

    def __init__(self, **kwargs):

        params = {}
        params['through_dir'] = 'LSST_THROUGHPUTS_BASELINE'
        params['atmos_dir'] = 'THROUGHPUTS_DIR'
        params['atmos'] = True
        params['aerosol'] = True
        params['telescope_files'] = ['detector.dat', 'lens1.dat',
                                     'lens2.dat', 'lens3.dat',
                                     'm1.dat', 'm2.dat', 'm3.dat']
        params['filterlist'] = 'ugrizy'
        params['wave_min'] = 300.
        params['wave_max'] = 1150.
        # params['load_components'] = False

        for par in ['through_dir', 'atmos_dir', 'atmos', 'aerosol',
                    'telescope_files', 'filterlist', 'wave_min',
                    'wave_max']:
            if par in kwargs.keys():
                params[par] = kwargs[par]
                # params[par]=str(kwargs[par])

        # self.atmos = params['atmos']
        # self.throughputsDir = os.getenv(params['through_dir'])
        self.throughputsDir = params['through_dir']

        self.atmosDir = params['atmos_dir']

        self.telescope_files = params['telescope_files']
        self.filter_files = ['filter_'+f+'.dat' for f in params['filterlist']]
        if 'filter_files' in kwargs.keys():
            self.filter_files = kwargs['filter_files']
        self.wave_min = params['wave_min']
        self.wave_max = params['wave_max']

        self.filterlist = params['filterlist']
        self.filtercolors = {'u': 'b', 'g': 'c',
                             'r': 'g', 'i': 'y', 'z': 'r', 'y': 'm'}

        self.lsst_std = {}
        self.lsst_system = {}
        self.mean_wavelength = {}
        self.lsst_detector = {}
        self.lsst_atmos = {}
        self.lsst_atmos_aerosol = {}
        # airmass = 1.0
        # aerosol = ''
        self.load_system()
        self.load_darkSky()
        self.load_components()

        """
        if params['atmos']:
            self.Load_Atmosphere(airmass=airmass, aerosol=aerosol)
        else:
            for f in self.filterlist:
                self.lsst_atmos[f] = self.lsst_system[f]
                self.lsst_atmos_aerosol[f] = self.lsst_system[f]

        """
        # self.lsst_telescope={}

        # self.Load_Telescope()

        # self.mean_wave()

    @property
    def system(self):
        return self.lsst_system

    @property
    def telescope(self):
        return self.lsst_telescope

    @property
    def atmosphere(self):
        return self.lsst_atmos

    @property
    def aerosol(self):
        return self.lsst_atmos_aerosol

    def load_components(self):
        """
        Method to load all components (detector, lenses, filters, mirrors)

        Returns
        -------
        None.

        """

        self.optical = {}

        for fi in self.telescope_files:
            vv = Bandpass()
            path_vv = os.path.join(self.throughputsDir, fi)
            vv.read_throughput(path_vv)
            nname = fi.split('.dat')[0]
            self.optical[nname] = vv

        self.filter = {}
        for f in self.filterlist:
            vv = Bandpass()
            path_vv = os.path.join(self.throughputsDir,
                                   'filter_{}.dat'.format(f))
            vv.read_throughput(path_vv)
            self.filter[f] = vv

    def load_system(self):
        """ Load files required to estimate throughputs
        """

        for f in self.filterlist:
            self.lsst_std[f] = Bandpass()
            self.lsst_system[f] = Bandpass()

            telfiles = ''
            if len(self.telescope_files) > 0:
                index = [i for i, x in enumerate(
                    self.filter_files) if f+'.dat' in x]
                telfiles = self.telescope_files+[self.filter_files[index[0]]]
            else:
                telfiles = self.filter_files

            self.lsst_system[f].read_throughput_list(telfiles,
                                                     root_dir=self.throughputsDir,
                                                     wavelen_min=self.wave_min,
                                                     wavelen_max=self.wave_max)

    """
    def Load_Telescope(self):
        for syst.lsst_telescope[system] = Bandpass()
            self.lsst_telescope[system].readThroughputList([system],
                                                        rootDir=self.throughputsDir,
                                                        wavelen_min=self.wave_min,
                                                        wavelen_max=self.wave_max)
     """

    def load_darkSky(self):
        """ Load DarkSky
        """
        self.darksky = Sed()
        self.darksky.read_sed_flambda(os.path.join(
            self.throughputsDir, 'darksky.dat'))

    def load_Atmosphere_deprecated(self, airmass=1.2, aerosol=True):
        """ Load atmosphere files
        and convolve with transmissions

        Parameters
        --------------
        airmass : float,opt
          airmass value
          Default : 1.2
        """
        self.airmass = airmass
        self.aerosol = aerosol
        if self.airmass > 0.:
            atmosphere = Bandpass()
            path_atmos = os.path.join(
                self.atmosDir, 'atmos_%d.dat' % (self.airmass*10))
            if not os.path.exists(path_atmos):
                path_atmos = 'atmos.dat'
            atmosphere.read_throughput(path_atmos)
            self.atmos = Bandpass(wavelen=atmosphere.wavelen, sb=atmosphere.sb)

            for f in self.filterlist:
                wavelen, sb = self.lsst_system[f].multiply_throughputs(
                    atmosphere.wavelen, atmosphere.sb)
                self.lsst_atmos[f] = Bandpass(wavelen=wavelen, sb=sb)

            if self.aerosol:
                atmosphere_aero = Bandpass()
                path_aero_atmos = os.path.join(
                    self.atmosDir, 'atmos_%d_aerosol.dat' % (self.airmass*10))
                atmosphere_aero.read_throughput(path_aero_atmos)
                self.atmos_aerosol = Bandpass(
                    wavelen=atmosphere_aero.wavelen, sb=atmosphere_aero.sb)

                for f in self.filterlist:
                    wavelen, sb = self.lsst_system[f].multiply_throughputs(
                        atmosphere_aero.wavelen, atmosphere_aero.sb)
                    self.lsst_atmos_aerosol[f] = Bandpass(
                        wavelen=wavelen, sb=sb)
        else:
            for f in self.filterlist:
                self.lsst_atmos[f] = self.lsst_system[f]
                self.lsst_atmos_aerosol[f] = self.lsst_system[f]

    def load_atmosphere(self, airmass=1.2, aerosol='aerosol'):
        """ Load atmosphere files
        and convolve with transmissions

        Parameters
        --------------
        airmass : float,opt
          airmass value
          Default : 1.2
        """
        self.airmass = airmass
        self.aerosol_v = aerosol
        self.aerosol_b = False

        if airmass > 0.:
            airmass_val = str(int(10*airmass))
            fName = ['atmos']+[airmass_val]
            # fName = 'atmos_{}.dat'.format(int(10*airmass))
            atmosphere = self.get_bandpass(fName)
            self.atmos = atmosphere
            # self.lsst_atmos = self.get_throughputs(atmosphere)

            if self.aerosol_v != '':
                self.aerosol_b = True
                # fName = 'atmos_{}_{}.dat'.format(int(10*airmass), aerosol)
                fName = ['atmos']+[airmass_val]+[aerosol]
                atmosphere_aerosol = self.get_bandpass(fName)
                self.atmos_aerosol = atmosphere_aerosol
                self.lsst_atmos_aerosol = self.get_throughputs(
                    atmosphere_aerosol)

    def get_bandpass(self, fName):
        """
        Method to grab the band pass corresponding to data in fName

        Parameters
        ----------
        fName : str
            file name to process.

        Returns
        -------
        atmos : Bandpass
            wavelength and throughputs.

        """

        atmosphere = Bandpass()
        fNameb = '{}.dat'.format('_'.join(fName))
        path_atmos = os.path.join(self.atmosDir, fNameb)
        if not os.path.exists(path_atmos):
            # print('Big problem: file', fNameb, 'not found')
            fName[1] = '25'
            fNameb = '{}.dat'.format('_'.join(fName))
            path_atmos = os.path.join(self.atmosDir, fNameb)
            # print('using', fNameb, 'instead')
        atmosphere.read_throughput(path_atmos)
        atmos = Bandpass(wavelen=atmosphere.wavelen, sb=atmosphere.sb)

        return atmos

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
        for f in self.filterlist:
            wavelen, sb = self.lsst_system[f].multiply_throughputs(
                bandpass.wavelen, bandpass.sb)
            through[f] = Bandpass(wavelen=wavelen, sb=sb)

        return through

    def plot_Throughputs(self, fig=None, ax=None):
        """ Plot the throughputs
        """
        # colors=['b','g','r','m','c',[0.8,0,0]]
        # style = [',', ',', ',', ',']

        if fig is None:
            fig, ax = plt.subplots(figsize=(15, 12))

        for i, band in enumerate(self.filterlist):

            ax.plot(self.lsst_system[band].wavelen,
                    self.lsst_system[band].sb,
                    linestyle='--', color=self.filtercolors[band],
                    label='%s - syst' % (band))
            """
            ax.plot(self.lsst_atmos[band].wavelen,
                    self.lsst_atmos[band].sb,
                    linestyle='-.', color=self.filtercolors[band],
                    label='%s - syst+atm' % (band))
            """
            if len(self.lsst_atmos_aerosol) > 0:
                ax.plot(self.lsst_atmos_aerosol[band].wavelen,
                        self.lsst_atmos_aerosol[band].sb,
                        linestyle='-', color=self.filtercolors[band],
                        label='%s - syst+atm+aero' % (band))

        ax.plot(self.atmos.wavelen, self.atmos.sb, color='k',
                label='X =%.1f atmos' % (self.airmass), linestyle='-')
        if len(self.lsst_atmos_aerosol) > 0:
            ax.plot(self.atmos_aerosol.wavelen, self.atmos_aerosol.sb,
                    color='k',
                    label='X =%.1f atm+aero' % (self.airmass),
                    linestyle='--')
        # plt.legend(loc=(0.85, 0.1), fontsize='smaller',
            # fancybox=True, numpoints=1)

        ax.legend(loc=(0.82, 0.1), fancybox=True, numpoints=1)

        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Sb (0-1)')
        ax.set_title('System throughput')
        ax.grid(visible=True)

    def plot_DarkSky(self, fig=None, ax=None):
        """ Plot darksky
        """
        # self.Load_DarkSky()

        if fig is None:
            fig, ax = plt.subplots()

        ax.plot(self.darksky.wavelen,
                self.darksky.flambda, 'k:', linestyle='-')
        ax.xlabel('Wavelength (nm)')
        ax.ylabel('flambda (ergs/cm$^2$/s/nm)')
        ax.title('Dark Sky SED')

    def plot_Throughputs_Spectrum(self, wavelength, fluxes, z, fig=None, ax=None):
        """ Plot throughput spectrum
        """
        if fig is None:
            fig, ax = plt.subplots()
        # style = [',', ',', ',', ',']

        for i, band in enumerate(self.filterlist):

            ax.plot(self.lsst_system[band].wavelen,
                    self.lsst_system[band].sb,
                    linestyle='-', color=self.filtercolors[band],
                    label='%s - system' % (band))

    def mean_wave(self):
        """ Estimate mean wave
        """
        for band in self.filterlist:
            self.mean_wavelength[band] = np.sum(
                self.lsst_atmos_aerosol[band].wavelen*self.lsst_atmos_aerosol[band].sb)\
                / np.sum(self.lsst_atmos_aerosol[band].sb)
