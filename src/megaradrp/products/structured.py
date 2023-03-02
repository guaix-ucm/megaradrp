#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Products of the Megara Pipeline"""


import numina.types.structured as structured
import numina.core.tagexpr as tagexpr

import megaradrp.datamodel
from megaradrp.datatype import MegaraDataType

class BaseStructuredCalibration(structured.BaseStructuredCalibration):
    DATATYPE = MegaraDataType.STRUCT_PROCESSED

    def __init__(self, instrument='unknown'):
        datamodel = megaradrp.datamodel.MegaraDataModel()
        super(BaseStructuredCalibration, self).__init__(instrument, datamodel)
        my_tag_table = self.datamodel.query_attrs
        objtags = [my_tag_table[t] for t in self.tag_names()]
        self.query_expr = tagexpr.query_expr_from_attr(objtags)
        self.names_t = self.query_expr.tags()
        self.names_f = self.query_expr.fields()

        self.total_fibers = 0
        self.missing_fibers = []
        self.error_fitting = []

    def __getstate__(self):
        st = super(BaseStructuredCalibration, self).__getstate__()

        keys = ["total_fibers", "missing_fibers", "error_fitting"]

        for key in keys:
            st[key] = self.__dict__[key]

        return st

    def validate(self, obj):
        """Validate objects with the TRACE_MAP schema"""
        import json
        from numina.util.jsonencoder import ExtEncoder
        import megaradrp.validators as valid

        super(BaseStructuredCalibration, self).validate(obj)

        checker = valid.check_as_datatype(self.DATATYPE)
        # We have to convert obj to a dictionary
        # FIXME: Dumping to a string and reloading. This can be done better
        res = json.dumps(obj.__getstate__(), cls=ExtEncoder)
        serialized = json.loads(res)
        return checker(serialized)