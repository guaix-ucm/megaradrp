

from numina.core import BaseRecipe

from ..loader import load_drp


def test_recipes_are_defined():

    current_drp = load_drp()

    assert 'default' in current_drp.pipelines

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            assert isinstance(recipe, BaseRecipe)


def test_recipes_have_tags():

    current_drp = load_drp()

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            qfields = recipe.tag_names()
            print('+', recipe)
            print('-', qfields)
    assert True