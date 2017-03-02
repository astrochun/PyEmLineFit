from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pyspeckit as psk

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
    infile0 = '/Users/cly/data/DEEP2/DR4/DEEP2_2D_Field1.f3.fits'

    data0, hdr0 = fits.getdata(infile0, header=True)

    x0 = hdr0['CRVAL1'] + hdr0['CDELT1']*np.arange(hdr0['NAXIS1'])
    y0 = data0[4]

    lines = np.array([4958.91, 5006.84])
    zspec = 0.7835

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
