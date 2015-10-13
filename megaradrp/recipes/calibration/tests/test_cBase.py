import inspect
from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe


def test_number_of_methods():
    lista = ['__call__', '__init__', 'buildRI', 'build_recipe_input', 'configure', 'create_input', 'create_result', 'hdu_creation', 'products', 'requirements', 'run', 'set_base_headers']
    methods = inspect.getmembers(MegaraBaseRecipe, predicate=inspect.ismethod)

    assert len(methods) is 12

    methods.sort()
    for key,value in methods:
        assert key in lista
        lista.pop(0)





if __name__ == '__main__':
    test_number_of_methods()