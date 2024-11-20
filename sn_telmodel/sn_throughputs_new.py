#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:50:53 2024

@author: philippe.gris@clermont.in2p3.fr
"""
from functools import wraps
from rubin_sim.phot_utils import Bandpass, Sed
from sn_telmodel.sn_telescope_new import Telescope
from sn_telmodel.sn_atmosphere import Atmos_Transmission
from rubin_sim.phot_utils import Bandpass

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
                 darksky_file='throughputs_1.9/baseline/darksky.dat'):
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

        self.load_atmosphere()
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
            fig, ax = plt.subplots(figsize=(15, 12))

        for i, band in enumerate(self.filter_list):

            ax.plot(self.tel_trans[band].wavelen,
                    self.tel_trans[band].sb,
                    linestyle='--', color=self.filter_colors[band],
                    label='%s - tel' % (band))
            """
            ax.plot(self.lsst_atmos[band].wavelen,
                    self.lsst_atmos[band].sb,
                    linestyle='-.', color=self.filtercolors[band],
                    label='%s - tel+atm' % (band))
            """
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
        """ Load DarkSky
        """
        self.darksky = Sed()
        self.darksky.read_sed_flambda(darksky_file)

    def plot_darksky(self, plt, fig=None, ax=None):
        """ Plot darksky
        """
        # self.Load_DarkSky()

        if fig is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        ax.plot(self.darksky.wavelen,
                self.darksky.flambda, 'k:', linestyle='-')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('flambda (ergs/cm$^2$/s/nm)')
        fig.suptitle('Dark Sky SED')
        ax.grid(visible=True)
