__author__ = 'Pica4x6'


def test_cDark():
    # dark_image
    from numina.core import init_drp_system, import_object
    from megaradrp.recipes.calibration.cDark import DarkRecipe
    from megaradrp.requirements import MasterBiasRequirement

    print "Test cDark"

    drps = init_drp_system()
    instrument = drps.get('MEGARA')
    pipeline = instrument.pipelines.get('default')
    recipe_fqn = pipeline.recipes.get('dark_image')
    RecipeClass = import_object(recipe_fqn)

    assert RecipeClass is DarkRecipe

    recipe = RecipeClass()
    biasrequirement = MasterBiasRequirement()

    recipereqs = RecipeClass.Requirements(biasrequirement)

    reciperesult = recipe.run(recipereqs)