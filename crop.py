# 
#  crop.py
#  Simulation Software
#  
#  Created by Alexander Rudy on 2011-11-08.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 


import AstroObject.AstroImage as AI

filename = "CCD-BB-st.fits"

def mask(self,left,top,right=None,bottom=None):
    """Masks the image by the distances provided"""
    if not right:
        right = left
    if not bottom:
        bottom = top
    shape  = self.states[self.statename].shape
    masked = self.states[self.statename].data[left:shape[0]-right,top:shape[1]-bottom]
    # LOG.debug("Masked masked and saved image")
    self.save(masked,"Masked")
    
AI.ImageObject.mask = mask

Image = AI.ImageObject()

Image.read(filename)

print Image.list()

Image.keep(filename)
Image.write(filename.rstrip(".fits")+"sm.fits")
# x,y = Image.frame().shape
# 
# x -= 2048
# y -= 2048
# 
# x /= 2
# y /= 2
# 
# print x,y
# 
# Image.mask(x,y)
# 
# Image.write(filename,clobber=True)