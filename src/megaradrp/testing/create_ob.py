from astropy.io import fits as fits
from numina.core import ObservationResult
from numina.instrument import assembly as asb
from numina.types.dataframe import DataFrame


def generate_bias(detector, number, temporary_path):
    from megaradrp.simulation.actions import simulate_bias
    from megaradrp.recipes.calibration.bias import BiasRecipe

    config_uuid = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
    date_obs = "2017-11-09T11:00:00.0"
    fs = [simulate_bias(detector) for i in range(number)]
    header = fits.Header()
    header["DATE-OBS"] = date_obs
    header["INSCONF"] = config_uuid
    header["INSTRUME"] = "MEGARA"
    header["VPH"] = "LR-U"
    header["INSMODE"] = "MOS"
    for aux in range(len(fs)):
        fits.writeto(
            f"{temporary_path}/bias_{aux}.fits", fs[aux], header=header, overwrite=True
        )

    fs = [f"{temporary_path}/bias_{i}.fits" for i in range(number)]

    ob = ObservationResult()
    ob.instrument = "MEGARA"
    ob.mode = "bias_image"

    pkg_paths = ["megaradrp.instrument.configs"]
    store = asb.load_paths_store(pkg_paths)
    insmodel = asb.assembly_instrument(store, config_uuid, date_obs, by_key="uuid")
    insmodel.configure_with_header(header)
    ob.configuration = insmodel
    ob.frames = [DataFrame(filename=f) for f in fs]

    recipe = BiasRecipe()
    ri = recipe.create_input(obresult=ob)
    return recipe.run(ri)
