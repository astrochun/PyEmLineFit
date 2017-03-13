from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pyspeckit as psk

from pyspeckit import models

infile0 = '/Users/cly/data/DEEP2/DR4/DEEP2_2D_Field1.f3.fits'

data0, hdr0 = fits.getdata(infile0, header=True)

x0 = hdr0['CRVAL1'] + hdr0['CDELT1']*np.arange(hdr0['NAXIS1'])
y0 = data0[4]

zspec = 0.7835


#class OII_gauss(xarr, r_peak, r_center, r_sigma, b_peak, b_sigma):
def OII_gauss(xarr, r_peak, r_center, r_sigma, b_peak, b_sigma):
    # + on 13/03/2017
    xarr = np.array(xarr)

    b_center = r_center * 1.00074

    r_gauss = r_peak * np.exp(-(xarr-r_center)**2/(2*r_sigma**2))
    b_gauss = b_peak * np.exp(-(xarr-b_center)**2/(2*b_sigma**2))

    #def __init__():
    #    self.npars = 4
    #    self.npeaks = 2

    return r_gauss + b_gauss

# + on 13/03/2017
OII_parnames = ['r_peak', 'r_center', 'r_sigma', 'b_peak', 'b_center']
parlimited=[(True,False),(True,True),(True,False),(True,False),(True,True)]
parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)]

# + on 13/03/2017
OII_gauss_fitter = models.model.SpectralModel(OII_gauss, 5, parnames=OII_parnames,
                                              parlimited=parlimited, parlimits=parlimits,
                                              fitunits='angstrom')


# OII_fitter = psk.models.model.SpectralModel(OIIFitter, 5,
def get_exclude(z_lines, w=5):
    # Mod on 02/03/2017 to improve efficiency for exclude list definition
    listoflist = [[a,b] for a,b in zip(z_lines-w,z_lines+w)]
    exclude = np.reshape(listoflist, len(z_lines)*2)
    #for zz in range(len(z_lines)):
    #    test = (z_lines[zz]+np.array([-5,5])).tolist()
    #    exclude.append(test)
    #exclude.append((z_lines-5.0).tolist())
    #exclude.append((z_lines+5.0).tolist())

    return exclude.tolist()

def test():
    lines = np.array([4958.91, 5006.84])

    z_lines = lines * (1+zspec)
    z_line  = z_lines[1]

    box_reg = np.where((x0 >= z_line - 100) & (x0 <= z_line + 100))[0]

    sig0 = np.std(y0[box_reg])
    # med0 = np.median(y0[box_reg])

    print sig0 #, med0
    px = psk.units.SpectroscopicAxis(x0[box_reg], unit='angstroms')

    y0 = y0 # - med0
    sp = psk.Spectrum(data=y0[box_reg], xarr=px, unit='erg/s/cm2/Ang',
                      error=np.repeat(sig0,len(box_reg)))

    print 'before : ', np.min(sp.data), np.max(sp.data)
    exclude = get_exclude(z_lines)
    print len(exclude)
    print exclude

    sp.baseline(annotate=True, subtract=False, exclude=exclude)
    print 'after : ', np.min(sp.data), np.max(sp.data)

    med0 = sp.baseline.basespec
    # print '## med0 : ', med0

    guess = [y0[box_reg].max(), z_line, 1.0]
    print guess
    sp.specfit(fittype='gaussian', guesses=guess)
    print sp.specfit.parinfo

    # print type(sp.specfit.parinfo)

    fig, ax = plt.subplots()
    print type(ax)
    sp.plotter(ax)

    sp.specfit.plot_fit()
    #print type(ax)
    #ax.axhline(med0, color='m', linewidth=4)

    ax.plot(x0,y0)

    sp.baseline.plot_baseline(annotate=True, linewidth=2)
    cont_arr = sp.baseline.get_model(px) #x0[box_reg])
    #ax.plot(x0[box_reg],cont_arr, 'g-', linewidth=10)

    ax.set_xlim([8500,9000])

    fig = plt.gcf()
    fig.savefig('test.pdf')
    #fig.close()

    #plt.show()

    return sp
#endef

def test_OII():
    lines = np.array([3726.16, 3728.91])

    z_lines = lines * (1+zspec)
    z_line  = z_lines[1]

    box_reg = np.where((x0 >= z_line - 100) & (x0 <= z_line + 100))[0]

    sig0 = np.std(y0[box_reg])

    print sig0 #, med0
    px = psk.units.SpectroscopicAxis(x0[box_reg], unit='angstroms')

    y0 = y0 # - med0
    sp = psk.Spectrum(data=y0[box_reg], xarr=px, unit='erg/s/cm2/Ang',
                      error=np.repeat(sig0,len(box_reg)))

    print 'before : ', np.min(sp.data), np.max(sp.data)
    exclude = get_exclude(z_lines)
    print len(exclude)
    print exclude

    sp.baseline(annotate=True, subtract=False, exclude=exclude)
    print 'after : ', np.min(sp.data), np.max(sp.data)

    med0 = sp.baseline.basespec

    sp.Registry.add_fitter('OII_gauss', OII_gauss_fitter, 5)
    sp.specfit.register_fitter(name='OII_gauss', function=OII_gauss_fitter, npars=5)

    guess = [np.max(y0[box_reg]), z_lines[1], 1.0, 0.7*np.max(y0[box_reg]), 1.0]
    #guess = [y0[box_reg].max(), z_line, 1.0]
    #print guess
    limits = [(0,0), (z_lines[1]-2,z_lines[1]+2), (0,5), (0,0), (0,5)]
    sp.specfit(fittype='OII_gauss', guesses=guess, limits=limits)
    print sp.specfit.parinfo

    fig, ax = plt.subplots()
    print type(ax)
    sp.plotter(ax)

    sp.specfit.plot_fit()

    ax.plot(x0,y0)

    sp.baseline.plot_baseline(annotate=True, linewidth=2)
    cont_arr = sp.baseline.get_model(px) #x0[box_reg])

    ax.set_xlim([6600,6700])

    fig = plt.gcf()
    fig.savefig('test_OII.pdf')
    #fig.close()

    #plt.show()

    return sp
#endef

def test_single():
    lines   = np.array([4861.32, 4958.91, 5006.84])
    z_lines = lines * (1+zspec)

    OIII5007 = z_lines[2]
    OIII4959 = z_lines[1]
    Hbeta    = z_lines[0]
    sig0 = np.std(y0)
    px = psk.units.SpectroscopicAxis(x0, unit='angstroms')

    sp = psk.Spectrum(data=y0, xarr=px, unit='erg/s/cm2/Ang',
                      error=np.repeat(sig0, len(y0)))

    exclude = get_exclude(z_lines)
    print len(exclude)
    sp.baseline(annotate=True, subtract=False, exclude=exclude)

    #Hbeta sp.specfit.selectregion(xmin=Hbeta-10, xmax=Hbeta+10)
    ampHb = np.max(sp.data[sp.specfit.xmin:sp.specfit.xmax])

    #OIII4959
    sp.specfit.selectregion(xmin=OIII4959-10, xmax=OIII4959+20)
    amp4959 = np.max(sp.data[sp.specfit.xmin:sp.specfit.xmax])

    #OIII5007
    sp.specfit.selectregion(xmin=OIII5007-10, xmax=OIII5007+20)
    amp5007 = np.max(sp.data[sp.specfit.xmin:sp.specfit.xmax])

    guess = [ampHb, Hbeta, 1.0, amp4959, OIII4959, 1.0,
             amp5007, OIII5007, 1.0]

    fig, ax = plt.subplots()
    sp.plotter(ax)

    sp.specfit(guesses=guess, negamp=False, quiet=True,
               show_components=True)
    print sp.specfit.parinfo

    sp.baseline.plot_baseline(annotate=True, linewidth=2)

    sp.plotter.axis.set_xlim(8600,9200)
    # sp.specfit.plot_fit()

    fig = plt.gcf()
    fig.savefig('test_single.pdf')

#enddef
