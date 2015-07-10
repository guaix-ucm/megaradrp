__author__ = 'Pica4x6'

def test_cBias():
    from numina.core import init_drp_system, import_object
    from numina.core import ObservationResult, DataFrame
    from megaradrp.recipes.calibration.cBias import BiasRecipe

    print "Test cBias"

    drps = init_drp_system()
    instrument = drps.get('MEGARA')
    pipeline = instrument.pipelines.get('default')
    recipe_fqn = pipeline.recipes.get('bias_image')
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is BiasRecipe

    recipe = RecipeClass()

    obsresult = ObservationResult(mode='bias_image')

    obsresult.mode='bias_image'
    obsresult.images=[DataFrame(filename='/home/pica/Documents/Megara/_data/r00002.fits')]
    obsresult.instrument = 'MEGARA'


    recipereqs = RecipeClass.Requirements(obresult=obsresult)

    reciperesult = recipe.run(recipereqs)