#!/usr/bin/env python
# 
#  test.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


import SEDModel.Spectra as Spectra
import matplotlib.pyplot as plt

import logging
LOG = logging.getLogger("Test Program")
LOG.info("Starting SEDMachine Model Test Program.")
LOG.debug("Done Importing Modules.")

try:
    sp = Spectra.Spectrum()
    sp.generate_ramp(370,1000,10)
    sp.show_spectrum()
    plt.show()
except:
    Spectra.LOG.exception("Critical Failure!")
    raise




LOG.info("Ending SEDMachine Model Test Program.")