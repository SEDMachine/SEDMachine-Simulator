# 
#  Test_SEDImage.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-01.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


from SEDImage import *

import AstroObject.Tests.Test_AstroObjectAPI as BaseAPITests

class Test_SEDImage(BaseAPITests.API_Base_Object):
    """SEDImage.SEDImage"""
    def setUp(self):
        """Set up the API Tests"""
        self.VALID = np.zeros((100,100))
        self.INVALID = "STRING"
        self.SHOWTYPE = mpl.image.AxesImage
        self.HDUTYPE = pf.ImageHDU
        self.imHDU = pf.ImageHDU(self.VALID)
        self.HDU = pf.PrimaryHDU(self.VALID)
        self.OBJECT = SEDImage
        self.FRAMELABEL = "Picture"
        self.FRAME = ImageFrame
        self.FRAMEINST = ImageFrame(self.VALID,self.FRAMELABEL)
        self.OBJECTSTR = None
        
        def SAME(first,second):
            """Return whether these two are the same"""
            return not (np.abs(first-second) > 1e-6).any()
            
        self.SAME = SAME
        
        self.check_constants()
        