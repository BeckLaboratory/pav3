"""Tests for filter and merge configuration defaults."""
from __future__ import annotations

import pytest

from pav3 import const
from pav3.align.lcmodel import get_model
from pav3.align.lcmodel._lcmodel_null import LCAlignModelNull
from pav3.params import PavParams
from pav3.workflow.call import get_merge_params


def test_filter_defaults():
    """LC filtering is off and MAPQ filtering is on by default."""
    params = PavParams('asm', {})

    assert params.lc_model == 'none'
    assert isinstance(get_model(params.lc_model), LCAlignModelNull)
    assert params.min_mapq == 20


def test_lc_model_named():
    """The packaged default model is still available by name."""
    assert not isinstance(get_model('default'), LCAlignModelNull)
    assert isinstance(get_model('null'), LCAlignModelNull)


def test_lc_model_unknown_raises():
    """An unknown LC model name is an error, not a silent null model."""
    with pytest.raises(ValueError, match='not found'):
        get_model('no_such_lc_model_name')


def test_min_mapq_config():
    assert PavParams('asm', {'min_mapq': '0'}).min_mapq == 0

    with pytest.raises(ValueError):
        _ = PavParams('asm', {'min_mapq': -1}).min_mapq


def test_merge_params_default():
    """Without an override, merge parameters are the defaults (copies, so callers cannot modify the defaults)."""
    params = get_merge_params('insdel')

    assert params == const.DEFAULT_MERGE_PARAMS['insdel']
    assert params[0] is not const.DEFAULT_MERGE_PARAMS['insdel'][0]
    assert get_merge_params('ins') == get_merge_params('del') == params


def test_merge_params_override():
    """An override replaces the stages for listed variant types only."""
    override = {'insdel': {'ro_min': 0.5}, 'snv': [{'offset_max': 0, 'match_alt': True}]}

    assert get_merge_params('insdel', override) == [{'ro_min': 0.5}]
    assert get_merge_params('snv', override) == [{'offset_max': 0, 'match_alt': True}]
    assert get_merge_params('cpx', override) == const.DEFAULT_MERGE_PARAMS['cpx']


def test_merge_params_config_json():
    """The config parameter accepts a JSON object or a JSON string (assembly table config column)."""
    assert PavParams('asm', {'merge_params': {'inv': {'ro_min': 0.5}}}).merge_params == {'inv': {'ro_min': 0.5}}
    assert PavParams('asm', {'merge_params': '{"inv": {"ro_min": 0.5}}'}).merge_params == {'inv': {'ro_min': 0.5}}
    assert PavParams('asm', {}).merge_params is None

    with pytest.raises(ValueError):
        _ = PavParams('asm', {'merge_params': '[1, 2]'}).merge_params


@pytest.mark.parametrize('override', [{'indel': [{'ro_min': 0.5}]}, {'insdel': []}, {'insdel': [0.5]}])
def test_merge_params_invalid(override):
    with pytest.raises(ValueError):
        get_merge_params('insdel', override)
