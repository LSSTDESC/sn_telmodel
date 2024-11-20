#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 18:50:53 2024

@author: philippe.gris@clermont.in2p3.fr
"""

from sn_telmodel.sn_telescope_new import Telescope
from sn_telmodel.sn_atmosphere import Atmos_Transmission


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
