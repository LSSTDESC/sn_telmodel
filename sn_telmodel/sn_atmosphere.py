#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:47:41 2024

@author: philippe.gris@clermont.in2p3.fr
"""
from getObsAtmo.getObsAtmo import ObsAtmo
from sn_telmodel.sn_transtools import get_trans
from rubin_sim.phot_utils import Bandpass
import numpy as np
import glob


class Atmos_Transmission:
    def __init__(self, site_name='LSST', pressure=743.,
                 atmos_dir='throughputs_1.9/atmos', atmos_type='obsatmo'):
        """
        class to estimate atmospheric transmission

        Parameters
        ----------
        site_name : str, optional
            Site name. The default is 'LSST'.
        pressure : float, optional
            Pressure. The default is 743..
        atmos_dir : str, optional
            Location dir of atmosphere files. 
            The default is 'throughput_v1.9/atmos'.
        atmos_type : str, optional
            Transmission estimation method . The default is 'obsatmo'.

        Returns
        -------
        None.

        """

        self.site_name = site_name
        self.pressure = pressure
        self.atmos_dir = atmos_dir
        self.atmos_type = atmos_type

        if self.atmos_type == 'obsatmo':
            self.emul = ObsAtmo(site_name, pressure)

    def load_atmosphere(self, site_name='LSST', airmass=1.2, aerosol=0.0,
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

        if self.atmos_type == 'obsatmo':
            self.load_atmosphere_obsatmo(site_name, airmass, aerosol,
                                         pwv, oz, beta, pressure)

        if self.atmos_type == 'from_file':
            self.load_atmosphere_from_file(airmass)

    def load_atmosphere_obsatmo(self, site_name='LSST', airmass=1.2, aerosol=0.0,
                                pwv=4.0, oz=300, beta=1.4, pressure=743.):
        """
        Load atmosphere files
        and convolve with transmissions

        Parameters
        ----------
        site_name : str, optional
            Site Name. The default is 'LSST'.
        airmass : float, optional
            airmass value. The default is 1.2.
        aerosol : float, optional
            aerosol value. The default is 0.0.
        pwv : float, optional
            precipitable water vapor value. The default is 4.0.
        oz : float, optional
            ozone value. The default is 300.
        beta : float, optional
            angström parameter. The default is 1.4.
        pressure : float, optional
            Pressure value. The default is 743..

        Returns
        -------
        None.

        """
        # emulate LSST
        if site_name != self.site_name or pressure != self.pressure:
            self.emul = ObsAtmo(site_name, pressure)
        self.atmosphere = self.get_bandpass_obsatmo(self.emul,
                                                    airmass, aerosol,
                                                    pwv, oz, beta)
        # self.atmos_aerosol = atmosphere_aerosol
        # self.lsst_atmos_aerosol = self.get_throughputs(atmosphere_aerosol)
        self.airmass = airmass
        self.aerosol = aerosol
        self.pwv = pwv
        self.oz = oz
        self.beta = beta
        self.pressure = pressure

    def get_bandpass_obsatmo(self, emul, airmass=1.2, aerosol=0.0,
                             pwv=4.0, oz=300, beta=1.4):
        """
        Method to get band pass using getObsAtmo

        Parameters
        ----------
        emul: getObsAtmo emulator
             used to estimate transmission
        airmass : float, optional
            Airmass value. The default is 1.2.
        aerosol : float, optional
            aerosol. The default is 0.0.
        pwv : float, optional
            precipitable water vapor. The default is 4.0.
        oz : float, optional
            Ozone. The default is 300.

        Returns
        -------
        atmos : Bandpass
            Atmospheric transmission.

        """

        trans = get_trans(airmass, pwv, oz, aerosol, beta,
                          colname=['wl', 'trans'],
                          emul=emul)
        trans = trans.round({'wl': 1, 'trans': 8})
        atmos = Bandpass(
            wavelen=np.asarray(trans['wl'].to_list()),
            sb=np.asarray(trans['trans'].to_list()))

        return atmos

    def load_atmosphere_from_file(self, airmass=1.2):
        """ Load atmosphere files
        and convolve with transmissions

        Parameters
        --------------
        airmass : float,opt
          airmass value
          Default : 1.2
        """

        fName = 'atmos_{}_aerosol'.format(int(10*airmass))

        fName_full = '{}/{}*.dat'.format(self.atmos_dir, fName)

        fis = glob.glob(fName_full)

        self.atmosphere = self.get_bandpass_from_file(fis[0])
        # self.lsst_atmos_aerosol = self.get_throughputs(atmosphere_aerosol)
        self.airmass = airmass
        self.aerosol = 0.04
        self.pwv = 4.0
        self.oz = 300.
        self.beta = 1.4
        self.pressure = 743.

    def get_bandpass_from_file(self, fName):
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
        atmosphere.read_throughput(fName)
        atmos = Bandpass(wavelen=atmosphere.wavelen, sb=atmosphere.sb)

        return atmos

    def plot_atmospheric_transmission(self, plt, fig=None, ax=None,
                                      figtit='Atmospheric transmission',
                                      label='', linestyle='solid', color='b'):
        """
        Method to plot atmospheric transmission

        Parameters
        ----------
        plt : matplotlib.pyplot
            To make the plot.
        fig : matplotlib figure, optional
            Figure for the plot. The default is None.
        ax : matplotlib axis, optional
            Axis for the plot. The default is None.

        Returns
        -------
        None.

        """

        if fig is None:
            fig, ax = plt.subplots(figsize=(12, 8))

        ax.plot(self.atmosphere.wavelen, self.atmosphere.sb, label=label,
                color=color, linestyle=linestyle)

        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Sb (0-1)')
        if figtit != '':
            ax.set_title(figtit)
        ax.set_ylim([0.0, 1.])
        ax.grid(visible=True)
        if label != '':
            ax.legend()
