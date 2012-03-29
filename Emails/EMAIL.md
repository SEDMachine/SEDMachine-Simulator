# Choong: Blackbody

Attached are the codes from today.
They rely on the following Python libraries: matplotlib (plotting), scipy
(interpolation) and numpy (arrays), all of which are important scientific
python libraries.
-Alex

# Nick: Simulation Code

Hi Nick,

I wanted to clarify some of the code I have from you, to make sure I am understanding all of the data correctly.

In go36.py, there are a few variables that I don't think I understand:

widthmm (this looks like the width of the focal plane in mm, as it is used to recenter the lenslet positions)

p1,p2 (columns 2 and 3 from the data file) appear to be pupil positions? In mm?

# Nick: Wavelength Dependent PSF

Hi Nick,

Sorry it has taken me so long to get this to you. I've been having some memory/crashing problems this week with the simulator, so I've done a little bit of work stabilizing the current status. I've got a plan going forwards for further stabilization and optimization, so I think I'm in a good position to conitnue developing new features. With that in mind, PSFs!

 Here are my thoughts on wavelength dependent PSFs:

1. Continuous PSFs have to be descritized at some level. I proposed discritizing them at the sub-pixel level (i.e. 5 times per larger pixel) to result in a roughly continuous PSF across the spectrum.
2. The PSFs will be pre-computed at a variety of wavelengths, and then the closest one will be applied. I'll make the number of PSFs calculated configurable, but I'd like to avoid calculating PSFs dynamically, but I'll leave that ability in the code for *really* accurate simulations.
3. If I have the PSF in encircled energy form, calculating the full PSF mask is easy. As such, any sort of table that has wavelength along one axis, and distance from the center of the PSF along the other, should do for the format.
4. Along with this though, there will have to be some interpolation between wavelengths to find the PSFs. Even if we specify the PSF very densly, we can't possibly specify it exactly for each spectrum. I'll need some sort of 2-dimensional interpolation to turn our sparse grid of PSFs into a function of some kind.

Here is the status of other things I'm working on:

1. Unit/Flux calibration of everything. I have the flux and spectra calculations included in the simulation code right now, but I have put no effort so far into ensuring that the units come out correct in the end. I'll try to work on that this week, it should simply be a matter of following the units though each stage and ensuring that I can provide configurable conversions at each step.
2. Dark and Bias calculations: Currently, there is some sense of dark and bias current, but again, the units are not very robust. The basic principle is that each one is using a random Poisson distribution of counts around the specified mean count.
3. Re-use of all aspects of calculation: I'm slowly working on ways to Cache more and more data between runs. This will make simulation across a variety of Temperatures etc. much faster, as the program will try to only re-calculate data as needed.

This stuff might require us to skype. If so, let me know.

-Alex
