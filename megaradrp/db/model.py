#
# Copyright 2017 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

"""Custom tables for MEGARA"""


from sqlalchemy import Integer, String, DateTime, Float
from sqlalchemy import CHAR
from sqlalchemy import Table, Column, ForeignKey
from sqlalchemy.orm import relationship, synonym


from numinadb.base import Base


class MegaraFrame(Base):
    __tablename__ = 'megara_frames'
    id = Column(Integer, primary_key=True)
    uuid = Column(CHAR(32), nullable=True)
    name = Column(String(100), nullable=False)
    # name = Column(String(100), unique=True, nullable=False)
    ob_id = Column(String,  ForeignKey("obs.id"), nullable=False)
    object = Column(String)
    start_time = Column(DateTime)
    exposure_time = Column(Float)
    insmode = Column(String)
    vph = Column(String)
    completion_time = Column(DateTime)
    # ob = relationship("MyOb", back_populates='frames')
    #
    # filename = synonym("name")


import datetime


def function(session, some, meta):
    newframe = MegaraFrame()
    print(some, meta)

    newframe.name = meta['path']
    newframe.uuid = meta['uuid']
    newframe.start_time = meta['observation_date']
    newframe.ob_id = meta['blckuuid']
    # No way of knowing when the readout ends...
    newframe.completion_time = newframe.start_time + datetime.timedelta(seconds=meta['darktime'])
    newframe.exposure_time = meta['exptime']
    newframe.object = meta['object']
    newframe.insmode = meta['insmode']
    newframe.vph = meta['vph']
    # ob.frames.append(newframe)
    # ob.object = meta['object']
    session.add(newframe)
