

def create_simple_hdul():
    """Create a simple image for testing"""
    import astropy.io.fits as fits

    prim = fits.PrimaryHDU()
    prim.header['instrume'] = 'MEGARA'
    prim.header['VPH'] = 'LR-B'
    prim.header['DATE-OBS'] = '2019-02-21T01:02:02.2'
    prim.header['insmode'] = 'LCB'
    prim.header['UUID'] = '410f6ea0-c3df-481c-9820-5cf9e0ed3d91'
    fibers = fits.ImageHDU(name='FIBERS')
    fibers.header['CONFID'] = 'a908e9d2-7599-41f3-9e9d-91042df01da5'
    simple_img = fits.HDUList([prim, fibers])
    return simple_img