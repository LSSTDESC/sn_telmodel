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

    return df
