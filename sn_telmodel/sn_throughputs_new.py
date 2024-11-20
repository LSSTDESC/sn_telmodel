#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:50:53 2024

@author: philippe.gris@clermont.in2p3.fr
"""

from sn_telmodel.sn_telescope_new import Telescope
from sn_telmodel.sn_atmosphere import Atmos_Transmission
from rubin_sim.phot_utils import Bandpass


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
                 atmos_type='obsatmo'):
        Telescope.__init__(self, tel_dir, tel_optical_files,
                           tel_filter_files, tel_wave_min, tel_wave_max,
                           filter_colors)
        Atmos_Transmission.__init__(self, site_name, pressure, atmos_dir,
                                    atmos_type)

        print('there man in Throughputs')
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
        """ Plot the throughputs
        """
        # colors=['b','g','r','m','c',[0.8,0,0]]
        # style = [',', ',', ',', ',']

        if fig is None:
            fig, ax = plt.subplots(figsize=(15, 12))

        for i, band in enumerate(self.filter_list):

            ax.plot(self.tel_trans[band].wavelen,
                    self.tel_trans[band].sb,
                    linestyle='--', color=self.filter_colors[band],
                    label='%s - syst' % (band))
            """
            ax.plot(self.lsst_atmos[band].wavelen,
                    self.lsst_atmos[band].sb,
                    linestyle='-.', color=self.filtercolors[band],
                    label='%s - syst+atm' % (band))
            """
            ax.plot(self.throughputs[band].wavelen,
                    self.throughputs[band].sb,
                    linestyle='-',
                    color=self.filter_colors[band],
                    label='%s - syst+atm+aero' % (band))

        # ax.plot(self.atmos.wavelen, self.atmos.sb, color='k',
        #        label='X =%.1f atmos' % (self.airmass), linestyle='-')

        ax.plot(self.atmosphere.wavelen, self.atmosphere.sb,
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
