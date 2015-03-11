#
# Copyright 2015 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#


from astropy.io import fits

def get_tags_from_full_ob(ob, reqtags=None):
    # each instrument should have one
    # perhaps each mode...
    files = ob.files
    cfiles = ob.children
    alltags = {}

    if reqtags is None:
        reqtags = []

    # Init alltags...
    # Open first image
    if files:
        for fname in files[:1]:
            with fits.open(fname) as fd:
                header = fd[0].header
                for t in reqtags:
                    alltags[t] = header[t]
    else:

        for prod in cfiles[:1]:
            prodtags = prod.tags
            for t in reqtags:
                alltags[t] = prodtags[t]
    
    for fname in files:
        with fits.open(fname) as fd:
            header = fd[0].header
            
            for t in reqtags:
                if alltags[t] != header[t]:
                    msg = 'wrong tag %s in file %s' % (thisfile[t], fname)
                    raise ValueError(msg)

    for prod in cfiles:
        prodtags = prod.tags
        for t in reqtags:
            if alltags[t] != prodtags[t]:
                msg = 'wrong tag %s in product %s' % (thisfile[t], prod)
                raise ValueError(msg)

    return alltags


def tagger_empty(obsres):
    return get_tags_from_full_ob(obsres, reqtags=[])

def tagger_vph(obsres):
    return get_tags_from_full_ob(obsres, reqtags=['vph'])

