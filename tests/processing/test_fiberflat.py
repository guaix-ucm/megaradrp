from astropy.io import fits
import numpy
import pytest

from megaradrp.datamodel import MegaraDataModel
from megaradrp.processing.fiberflat import CommonFlatCorrector, FiberFlatCorrector, FlipLR, Splitter


def make_img(value=10.0, uuid="410f6ea0-c3df-481c-9820-5cf9e0ed3d91"):
    hdu = fits.PrimaryHDU(numpy.full((5, 8), value, dtype="float32"))
    hdu.header["UUID"] = uuid
    return fits.HDUList([hdu])


def make_flat():
    flat = numpy.full((5, 8), 2.0, dtype="float32")
    flat[1, 3] = 0.0
    return flat


@pytest.mark.parametrize("kind", ["array", "hdu", "hdulist"])
def test_common_flat_corrector(kind):
    flat = make_flat()
    if kind == "hdu":
        flatdata = fits.ImageHDU(flat)
    elif kind == "hdulist":
        flatdata = fits.HDUList([fits.PrimaryHDU(flat)])
    else:
        flatdata = flat

    corrector = CommonFlatCorrector(flatdata, calibid="flat-id")
    img = corrector(make_img())

    expected = numpy.full((5, 8), 5.0)
    # the pixels where the flat is zero are not corrected
    expected[1, 3] = 10.0
    assert numpy.allclose(img[0].data, expected)
    assert img[0].header["NUM-FLT"] == "flat-id"
    assert "Flat correction 410f6ea0-c3df-481c-9820-5cf9e0ed3d91" in img[0].header["HISTORY"]


def test_fiber_flat_corrector():
    flat = fits.HDUList([fits.PrimaryHDU(make_flat())])
    flat[0].header["UUID"] = "11111111-2222-3333-4444-555555555555"
    corrector = FiberFlatCorrector(flat, datamodel=MegaraDataModel())
    img = corrector(make_img())

    assert numpy.allclose(img[0].data[0], 5.0)
    # the calibration is identified by the UUID of the fiber flat
    assert img[0].header["NUM-FIBF"] == "11111111-2222-3333-4444-555555555555"
    assert "Fiber flat correction 410f6ea0-c3df-481c-9820-5cf9e0ed3d91" in img[0].header["HISTORY"]


def test_fiber_flat_corrector_invalid_fibers():
    """The fibers invalid in the fiber flat are invalid in the corrected image"""
    from megaradrp.datamodel import create_default_fiber_header

    flat = fits.HDUList([fits.PrimaryHDU(make_flat()), fits.ImageHDU(name="FIBERS")])
    flat[0].header["UUID"] = "11111111-2222-3333-4444-555555555555"
    flat["FIBERS"].header = create_default_fiber_header("LCB")
    flat["FIBERS"].header["FIB003_V"] = False
    flat["FIBERS"].header["FIB003QF"] = "broken"
    img = make_img()
    img.append(fits.ImageHDU(name="FIBERS"))
    img["FIBERS"].header = create_default_fiber_header("LCB")
    img = FiberFlatCorrector(flat, datamodel=MegaraDataModel())(img)
    assert img["FIBERS"].header["FIB003_V"] is False
    assert "FIB004_V" not in img["FIBERS"].header


def zero_order_flat():
    flat = fits.HDUList([fits.PrimaryHDU(make_flat())])
    flat[0].header["UUID"] = "11111111-2222-3333-4444-555555555555"
    mask = numpy.zeros((5, 8), dtype="uint8")
    mask[:, 2:4] = 1
    mask[4, 6] = 2
    flat.append(fits.ImageHDU(mask, name="ZOMASK"))
    return flat, mask


def test_fiber_flat_corrector_zero_order_mask():
    """The mask of the zero order of the fiber flat is added to the corrected image"""
    flat, mask = zero_order_flat()
    img = FiberFlatCorrector(flat, datamodel=MegaraDataModel())(make_img())
    assert "ZOMASK" in img
    assert numpy.array_equal(img["ZOMASK"].data, mask)
    # the mask of the flat is not shared with the image
    img["ZOMASK"].data[:] = 0
    assert flat["ZOMASK"].data.any()


def test_fiber_flat_corrector_zero_order_mask_combined():
    """An image that already has a zero order mask keeps it, combined with OR"""
    flat, mask = zero_order_flat()
    img = make_img()
    prev = numpy.zeros((5, 8), dtype="uint8")
    prev[0, 0] = 4
    prev[0, 2] = 2
    img.append(fits.ImageHDU(prev, name="ZOMASK"))
    img = FiberFlatCorrector(flat, datamodel=MegaraDataModel())(img)
    assert len([hdu for hdu in img if hdu.name == "ZOMASK"]) == 1
    assert img["ZOMASK"].data[0, 0] == 4
    assert img["ZOMASK"].data[0, 2] == 3
    assert img["ZOMASK"].data[4, 6] == 2


def test_fiber_flat_corrector_zero_order_mask_shape(caplog):
    """A mask of another shape is not added"""
    flat, _ = zero_order_flat()
    flat["ZOMASK"].data = numpy.zeros((5, 9), dtype="uint8")
    img = FiberFlatCorrector(flat, datamodel=MegaraDataModel())(make_img())
    assert "ZOMASK" not in img
    assert "does not match" in caplog.text


def test_splitter():
    splitter = Splitter()
    img = make_img()
    result = splitter(img)

    assert result is img
    assert splitter.out is not img
    # the copy does not change with the image
    img[0].data[:] = 0
    assert numpy.allclose(splitter.out[0].data, 10.0)


def test_flip_lr():
    img = make_img()
    img[0].data[:] = numpy.arange(8)
    result = FlipLR()(img)
    assert numpy.allclose(result[0].data[0], numpy.arange(8)[::-1])
