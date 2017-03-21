def running_std(x0, y0, line_flag, OH_flag0, window=100.0):
    sig0_arr = np.zeros(len(x0))

    for ss in xrange(len(x0)):
        box_reg = np.where((x0 >= x0[ss]-window) & (x0 <= x0[ss]+window) &
                           (line_flag == 0) & (OH_flag0 == 0) &
                           (y0 != 0.0))[0]
        sig0_arr[ss] = np.std(y0[box_reg])
    return sig0_arr
#enddef

sig0_arr = running_std(x0, y0, line_flag, OH_flag0, window=100.0)
sp   = psk.Spectrum(data=y0, xarr=px, error=sig0_arr)
