

from numina.core import BaseRecipe

from ..loader import load_drp


def test_recipes_are_defined():

    current_drp = load_drp()

    assert 'default' in current_drp.pipelines

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            assert isinstance(recipe, BaseRecipe)
