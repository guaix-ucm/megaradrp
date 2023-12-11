#
# Copyright 2017-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Custom tables for MEGARA"""


from sqlalchemy import Integer, String, DateTime, Float
from sqlalchemy import CHAR
from sqlalchemy import Column, ForeignKey


from numinadb.base import Base


class MegaraFrame(Base):
    __tablename__ = 'megara_frames'
    id = Column(Integer, primary_key=True)
    uuid = Column(CHAR(32), nullable=True)
    name = Column(String(100), nullable=False, unique=True)
    # name = Column(String(100), unique=True, nullable=False)
    ob_id = Column(String,  ForeignKey("obs.id"), nullable=False)
    object = Column(String)
    start_time = Column(DateTime)
    exposure_time = Column(Float)
    insmode = Column(String)
    vph = Column(String)
    completion_time = Column(DateTime)
