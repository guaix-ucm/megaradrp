from astropy.io import fits as fits

from megaradrp.testing.common import make_obresult, raw_header


def generate_bias(detector, number, temporary_path):
    from megaradrp.simulation.actions import simulate_bias
    from megaradrp.recipes.calibration.bias import BiasRecipe

    header = raw_header()
    names = [f"{temporary_path}/bias_{idx}.fits" for idx in range(number)]
    for name in names:
        fits.writeto(name, simulate_bias(detector), header=header, overwrite=True)

    ob = make_obresult("bias_image", names, header)
    recipe = BiasRecipe()
    ri = recipe.create_input(obresult=ob)
    return recipe.run(ri)
