from __future__ import division
from __future__ import print_function

import argparse
import numpy as np
from numpy.polynomial import Polynomial
import sys
import numina.instrument.assembly as asb
import numina.types.structured as structured
from numina.array.display.ximshow import ximshow_file
from numina.array.display.pause_debugplot import pause_debugplot


def assign_boxes_to_fibers(pseudo_slit_config, insmode):
    """Read boxes in configuration file and assign values to fibid

    Parameters
    ----------
    pseudo_slit_config : dict
        Contains the association of fibers and boxes
    insmode : string
        Value of the INSMODE keyword: 'LCB' or 'MOS'.

    Returns
    -------
    fibid_with_box : list of strings
        List with string label that contains both the fibid and the
        box name.

    """
    fibid_with_box = []
    n1 = 1
    list_to_print = []
    for dumbox in pseudo_slit_config:
        nfibers = dumbox['nfibers']
        name = dumbox['name']
        n2 = n1 + nfibers
        fibid_with_box += \
            ["{}  [{}]".format(val1, val2)
             for val1, val2 in zip(range(n1, n2), [name] * nfibers)]
        dumstr ='Box {:>2},  fibers {:3d} - {:3d}'.format(name, n1, n2 - 1)
        list_to_print.append(dumstr)
        n1 = n2
    print('\n* Fiber description for INSMODE={}'.format(insmode))
    for dumstr in reversed(list_to_print):
        print(dumstr)
    print('---------------------------------')

    return fibid_with_box


def plot_aper(ax, center_model, xmin, xmax, ix_offset,
               rawimage, fibids, fiblabel, colour, correction=0):
    if xmin == xmax == 0:
        num = 4096
        xp = np.linspace(start=1, stop=4096, num=num)
    else:
        num = int(float(xmax - xmin + 1) + 0.5)
        xp = np.linspace(start=xmin, stop=xmax, num=num)
    yp = center_model(xp) + correction
    if rawimage:
        lcut = (yp > 2056.5)
        yp[lcut] += 100
    ax.plot(xp + ix_offset, yp + 1, color=colour, linestyle='dotted')
    if fibids:
        if xmin == xmax == 0:
            xmidpoint = 2048
        else:
            xmidpoint = (xmin+xmax)/2
        ax.text(xmidpoint, yp[int(num / 2)], fiblabel, fontsize=6,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="grey", ),
                color=colour, fontweight='bold', backgroundcolor='white',
                ha='center')


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description="description: overplot traces"
    )
    # positional parameters
    parser.add_argument("fits_file",
                        help="FITS image containing the spectra",
                        type=argparse.FileType('r'))
    parser.add_argument("traces_file",
                        help="JSON file with fiber traces",
                        type=argparse.FileType('r'))
    # optional parameters
    parser.add_argument("--rawimage",
                        help="FITS file is a RAW image (otherwise trimmed "
                             "image is assumed)",
                        action="store_true")
    parser.add_argument("--global_offset", "--global-offset",
                        nargs='+', type=float,
                        help="Global offset polynomial coefficients "
                             "(+upwards, -downwards)")
    parser.add_argument("--fibids",
                        help="Display fiber identification number",
                        action="store_true")
    parser.add_argument("--verbose",
                        help="Enhance verbosity",
                        action="store_true")
    parser.add_argument("--z1z2",
                        help="tuple z1,z2, minmax or None (use zscale)")
    parser.add_argument("--bbox",
                        help="bounding box tuple: nc1,nc2,ns1,ns2")
    parser.add_argument("--keystitle",
                        help="tuple of FITS keywords.format: " +
                             "key1,key2,...keyn.'format'")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy",
                        default="0,0,640,480")
    parser.add_argument("--pdffile",
                        help="ouput PDF file name",
                        type=argparse.FileType('w'))
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args(args=args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

    # global_offset in command line
    if args.global_offset is None:
        args_global_offset_set = False
        args.global_offset = [0.0]
    else:
        args_global_offset_set = True

    args_global_offset = Polynomial(args.global_offset)

    # read pdffile
    if args.pdffile is not None:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(args.pdffile.name)
    else:
        pdf = None

    ax = ximshow_file(args.fits_file.name,
                      args_cbar_orientation='vertical',
                      args_z1z2=args.z1z2,
                      args_bbox=args.bbox,
                      args_keystitle=args.keystitle,
                      args_geometry=args.geometry,
                      pdf=pdf,
                      show=False)

    # trace offsets for RAW images
    if args.rawimage:
        ix_offset = 51
    else:
        ix_offset = 1

    # read and display traces from JSON file
    # TODO: some checks, this should be done by validating the struct
    with open(args.traces_file.name, mode='r') as fd:
        import json
        data = json.load(fd)
        if 'type_fqn' not in data:
            raise ValueError("malformed JSON file, 'type_fqn' missing")
    #
    apers = structured.open(args.traces_file.name)
    # Load metadata from the traces
    meta_info = apers.meta_info

    origin = meta_info['origin']
    insconf_uuid = origin['insconf_uuid']
    # None is allowed
    date_obs = origin.get('date_obs')

    tags = apers.tags
    insmode = tags['insmode']

    # create instrument model
    pkg_paths = ['megaradrp.instrument.configs']
    store = asb.load_paths_store(pkg_paths)

    insmodel = asb.assembly_instrument(store, insconf_uuid, date_obs, by_key='uuid')

    pseudo_slit_config = insmodel.get_value('pseudoslit.boxes', **tags)

    fibid_with_box = list(assign_boxes_to_fibers(pseudo_slit_config, insmode))
    total_fibers = apers.total_fibers
    if total_fibers != len(fibid_with_box):
        raise ValueError('Mismatch between number of fibers and '
                         'expected number from account from boxes')

    if args_global_offset_set:
        global_offset = args_global_offset
    else:
        global_offset = apers.global_offset

    print('>>> Using global_offset:', global_offset)
    ref_column = apers.ref_column

    for geot in apers.contents:
        fibid = geot.fibid
        fiblabel = fibid_with_box[fibid - 1]
        start = geot.start
        stop = geot.stop
        # skip fibers without trace
        if geot.valid:
            center_model = geot.aper_center()
            y_at_ref_column = center_model(ref_column)
            correction = global_offset(y_at_ref_column)
            plot_aper(ax, center_model, start, stop, ix_offset, args.rawimage,
                       args.fibids, fiblabel, colour='blue', correction=correction)
        else:
            print('Warning ---> Missing fiber:', fibid_with_box[fibid - 1])

    if pdf is not None:
        pdf.savefig()
        pdf.close()
    else:
        pause_debugplot(12, pltshow=True, tight_layout=True)


if __name__ == "__main__":

    main()
