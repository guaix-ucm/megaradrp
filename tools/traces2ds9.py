from __future__ import division
from __future__ import print_function

import argparse
import json

from numina.core import import_object
from numina.util.jsonencoder import ExtEncoder


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='traces2ds9')
    # positional parameters
    parser.add_argument("json_file",
                        help="JSON file with fiber traces",
                        type=argparse.FileType('r'))
    parser.add_argument("ds9_file",
                        help="Output region file in ds9 format",
                        type=argparse.FileType('w'))
    # optional parameters
    parser.add_argument("--numpix",
                        help="Number of pixels/trace (default 100)",
                        default=100, type=int)
    parser.add_argument("--yoffset",
                        help="Vertical offset (+upwards, -downwards)",
                        default=0, type=float)
    parser.add_argument("--new_json",
                        help="New JSON file after applying specified yoffset",
                        type=argparse.FileType('w'))
    parser.add_argument("--rawimage",
                        help="FITS file is a RAW image (RSS assumed instead)",
                        action="store_true")
    parser.add_argument("--fibid_at",
                        help="Display fiber identification number at location",
                        default=0, type=int)

    args = parser.parse_args(args=args)

    # Rebuild product
    data = json.load(args.json_file)
    type_fqn = data.get('type_fqn', 'megaradrp.products.tracemap.TraceMap')
    klass = import_object(type_fqn)
    obj = klass.__new__(klass)
    obj.__setstate__(data)

    obj.global_offset += args.yoffset

    obj.to_ds9_reg(args.ds9_file, rawimage=args.rawimage,
                   numpix=args.numpix, fibid_at=args.fibid_at)

    if args.new_json is not None:
        json.dump(obj.__getstate__(), args.new_json, indent=2, cls=ExtEncoder)



if __name__ == "__main__":

    main()
