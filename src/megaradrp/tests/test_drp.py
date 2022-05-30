
import pytest

from numina.core import BaseRecipe
import numina.core.tagexpr as tagexpr
import numina.types.multitype as mt

from ..loader import load_drp


@pytest.fixture
def current_drp():
    return load_drp()


def simple_tagger(depos, keys):
    result = {}
    for k in keys:
        result[k] = depos[k]
    return result


def test_recipes_are_defined(current_drp):

    assert 'default' in current_drp.pipelines

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            assert isinstance(recipe, BaseRecipe)

expected_tags = {
    "MegaraSuccess": set(), "MegaraFail": set(), "MegaraBadPixelMask": set(),
    "MegaraBiasImage": set(), "MegaraDarkImage": set(),
    "MegaraArcCalibration": {'insmode', 'speclamp', 'vph'},
    "MegaraSlitFlat": set(), "MegaraTraceMap": set(),
    "MegaraModelMap": {'insmode', 'vph'},
    "MegaraFiberFlatImage": {'insmode', 'vph'},
    "MegaraTwilightFlatImage": {'confid', 'insmode', 'vph'},
    "MegaraFocusSpectrograph": {'insmode', 'vph'},
    "MegaraFocusTelescope": {'confid', 'insmode', 'vph'},
    "MegaraLcbAcquisition": {'confid', 'insmode', 'vph'},
    "MegaraLcbImage": {'confid', 'insmode', 'vph'},
    "MegaraLcbStdStar": {'confid', 'insmode', 'vph'},
    "MegaraMosAcquisition": {'confid', 'insmode', 'vph'},
    "MegaraMosImage": {'confid', 'insmode', 'vph'},
    "MegaraMosStdStar": {'confid', 'insmode', 'vph'},
    "MegaraExtinctionStar": set(), "MegaraSensitivityStar": set()
}


def test_recipes_have_tags(current_drp):

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            qfields = recipe.tag_names()
            assert qfields == expected_tags[key]


results1 = {
    'MegaraSuccess': {},
    'MegaraFail': {},
    'MegaraBadPixelMask': {'master_bias': 105},
    'MegaraBiasImage': {'master_bpm': 205},
    'MegaraDarkImage': {'master_bias': 105},
    'MegaraArcCalibration': {'master_bias': 105, 'master_bpm': 205, 'master_apertures': 11, 'lines_catalog': 1027},
    'MegaraSlitFlat': {'master_bpm': 205, 'master_bias': 105},
    'MegaraTraceMap': {'master_bias': 105, 'master_bpm': 205},
    'MegaraModelMap': {'master_bpm': 205, 'master_bias': 105, 'master_slitflat': 1, 'master_traces': 11},
    'MegaraFiberFlatImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_apertures': 11, 'master_wlcalib': 21},
    'MegaraTwilightFlatImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_apertures': 11, 'master_wlcalib': 21, 'master_fiberflat': 49},
    'MegaraFocusSpectrograph': {'master_bias': 105, 'master_bpm': 205, 'master_apertures': 11, 'master_wlcalib': 21},
    'MegaraFocusTelescope': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraLcbAcquisition': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraLcbImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraLcbStdStar': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraMosAcquisition': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraMosImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraMosStdStar': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 1, 'master_wlcalib': 21, 'master_fiberflat': 49, 'master_twilight': 59, 'master_apertures': 11},
    'MegaraExtinctionStar': {},
    'MegaraSensitivityStar': {},
}


results2 = {
    'MegaraSuccess': {},
    'MegaraFail': {},
    'MegaraBadPixelMask': {'master_bias': 105},
    'MegaraBiasImage': {'master_bpm': 205},
    'MegaraDarkImage': {'master_bias': 105},
    'MegaraArcCalibration': {'master_bias': 105, 'master_bpm': 205, 'master_apertures': 12, 'lines_catalog': 1027},
    'MegaraSlitFlat': {'master_bpm': 205, 'master_bias': 105},
    'MegaraTraceMap': {'master_bias': 105, 'master_bpm': 205},
    'MegaraModelMap': {'master_bpm': 205, 'master_bias': 105, 'master_slitflat': 2, 'master_traces': 12},
    'MegaraFiberFlatImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_apertures': 12, 'master_wlcalib': 22},
    'MegaraTwilightFlatImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_apertures': 12, 'master_wlcalib': 22, 'master_fiberflat': 42},
    'MegaraFocusSpectrograph': {'master_bias': 105, 'master_bpm': 205, 'master_apertures': 12, 'master_wlcalib': 22},
    'MegaraFocusTelescope': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraLcbAcquisition': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraLcbImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraLcbStdStar': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraMosAcquisition': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraMosImage': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraMosStdStar': {'master_bias': 105, 'master_bpm': 205, 'master_slitflat': 2, 'master_wlcalib': 22, 'master_fiberflat': 42, 'master_twilight': 52, 'master_apertures': 12},
    'MegaraExtinctionStar': {},
    'MegaraSensitivityStar': {},
}


ob_repo1 = {'vph': 'HR-I', 'insmode': 'MOS', 'confid': '123', 'speclamp': 'ThNe'}
ob_repo2 = {'vph': 'HR-I', 'insmode': 'LCB', 'confid': '000', 'speclamp': 'ThNe'}


@pytest.mark.parametrize("ob_repo, results", [
    (ob_repo1, results1), (ob_repo2, results2)
])
def test_recipes_extract_tags(current_drp, ob_repo, results):

    calibs = {
        'lines_catalog': [
            {'id': 1025, 'tags': {'vph': 'LR-I', 'speclamp': 'ThNe'}},
            {'id': 1026, 'tags': {'vph': 'LR-I', 'speclamp': 'ThAr'}},
            {'id': 1027, 'tags': {'vph': 'HR-I', 'speclamp': 'ThNe'}},
            {'id': 1027, 'tags': {'vph': 'HR-I', 'speclamp': 'ThAr'}},
        ],
        'master_bias': [
            {'id': 105, 'tags': {}},
            {'id': 101, 'tags': {}},
            {'id': 102, 'tags': {}}
        ],
        'master_dark': [
        ],
        'master_bpm': [
            {'id': 205, 'tags': {}},
            {'id': 201, 'tags': {}},
            {'id': 202, 'tags': {}}
        ],
        'master_slitflat': [
            {'id': 5, 'tags': {'vph': 'LR-I', 'insmode': 'MOS'}},
            {'id': 1, 'tags': {'vph': 'HR-I', 'insmode': 'MOS'}},
            {'id': 2, 'tags': {'vph': 'HR-I', 'insmode': 'LCB'}}
        ],
        'master_traces': [
            {'id': 15, 'tags': {'vph': 'LR-I', 'insmode': 'MOS'}},
            {'id': 11, 'tags': {'vph': 'HR-I', 'insmode': 'MOS'}},
            {'id': 12, 'tags': {'vph': 'HR-I', 'insmode': 'LCB'}}
        ],
        'master_wlcalib': [
            {'id': 25, 'tags': {'vph': 'LR-I', 'insmode': 'MOS'}},
            {'id': 21, 'tags': {'vph': 'HR-I', 'insmode': 'MOS'}},
            {'id': 22, 'tags': {'vph': 'HR-I', 'insmode': 'LCB'}}
        ],
        'master_fiberflat': [
            {'id': 45, 'tags': {'vph': 'LR-I', 'insmode': 'MOS', 'confid': '123'}},
            {'id': 49, 'tags': {'vph': 'HR-I', 'insmode': 'MOS', 'confid': '123'}},
            {'id': 41, 'tags': {'vph': 'HR-I', 'insmode': 'MOS', 'confid': '321'}},
            {'id': 42, 'tags': {'vph': 'HR-I', 'insmode': 'LCB', 'confid': '000'}}
        ],
        'master_twilight': [
            {'id': 55, 'tags': {'vph': 'LR-I', 'insmode': 'MOS', 'confid': '123'}},
            {'id': 59, 'tags': {'vph': 'HR-I', 'insmode': 'MOS', 'confid': '123'}},
            {'id': 51, 'tags': {'vph': 'HR-I', 'insmode': 'MOS', 'confid': '321'}},
            {'id': 52, 'tags': {'vph': 'HR-I', 'insmode': 'LCB', 'confid': '000'}}
        ]
    }

    calibs['master_apertures']  = calibs['master_traces']

    for pipeval in current_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)

            qfields = recipe.tag_names()

            ob_tags = simple_tagger(ob_repo, qfields)
            match_per_recipe = {}
            for rkey, val in recipe.requirements().items():

                if val.type.isproduct():
                    if isinstance(val.type, mt.MultiType):
                        try_these = val.type.node_type
                    else:
                        try_these = [val.type]

                    for vtype in try_these:
                        expr2 = vtype.query_expr.fill_placeholders(**ob_tags)
                        # Load table of products
                        try:
                            table = calibs[rkey]
                            for entry in table:
                                match = expr2.eval(**entry['tags'])
                                if isinstance(match, tagexpr.Expression):
                                    msg = f'{rkey}  tags are not complete: {match}'
                                    raise ValueError(msg)

                                if match == True:
                                    match_per_recipe[rkey] = entry['id']
                                    break
                            else:
                                pass
                        except KeyError:
                            pass
            assert results[key] == match_per_recipe
