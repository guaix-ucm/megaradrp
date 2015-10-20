import inspect
import sys
from megaradrp.recipes.calibration.cBase import MegaraBaseRecipe


def test_number_of_methods():

    if sys.version_info>(3,0):
        lista = ['RecipeInput', 'RecipeResult', '__call__', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'buildRI', 'build_recipe_input', 'configure', 'create_input', 'create_result', 'hdu_creation', '_MegaraBaseRecipe__generate_flow', 'products', 'requirements', 'run', 'set_base_headers', 'bias_process_common']
    else:
        lista = ['RecipeInput', 'RecipeResult', '_MegaraBaseRecipe__generate_flow', '__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'bias_process_common', 'buildRI', 'build_recipe_input', 'configure', 'create_input', 'create_result', 'hdu_creation', 'products', 'requirements', 'run', 'set_base_headers']
    methods = inspect.getmembers(MegaraBaseRecipe)

    assert len(methods) is len(lista)

    lista.sort()
    methods.sort()
    for key,value in methods:
        assert key in lista
        lista.pop(0)


if __name__ == '__main__':
    test_number_of_methods()