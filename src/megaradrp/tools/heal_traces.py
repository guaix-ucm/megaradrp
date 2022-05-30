from __future__ import division
from __future__ import print_function

import argparse
from copy import deepcopy
import json
import numpy as np
from numpy.polynomial import Polynomial
import sys
from uuid import uuid4

import numina.instrument.assembly as asb
from numina.array.display.polfit_residuals import polfit_residuals
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
            [f"{val1}  [{val2}]"
             for val1, val2 in zip(range(n1, n2), [name] * nfibers)]
        dumstr =f'Box {name:>2},  fibers {n1:3d} - {n2 - 1:3d}'
        list_to_print.append(dumstr)
        n1 = n2
    print(f'\n* Fiber description for INSMODE={insmode}')
    for dumstr in reversed(list_to_print):
        print(dumstr)
    print('---------------------------------')

    return fibid_with_box


def plot_trace(ax, coeff, xmin, xmax, ix_offset,
               rawimage, fibids, fiblabel, colour):
    if xmin == xmax == 0:
        num = 4096
        xp = np.linspace(start=1, stop=4096, num=num)
    else:
        num = int(float(xmax - xmin + 1) + 0.5)
        xp = np.linspace(start=xmin, stop=xmax, num=num)
    ypol = Polynomial(coeff)
    yp = ypol(xp)
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
        description="description: heal traces"
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
    parser.add_argument("--global_offset",
                        help="Global offset polynomial coefficients "
                             "(+upwards, -downwards)")
    parser.add_argument("--fibids",
                        help="Display fiber identification number",
                        action="store_true")
    parser.add_argument("--verbose",
                        help="Enhance verbosity",
                        action="store_true")
    parser.add_argument("--healing",
                        help="JSON healing file to improve traces",
                        type=argparse.FileType('r'))
    parser.add_argument("--updated_traces",
                        help="JSON file with modified fiber traces",
                        type=argparse.FileType('w'))
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
        args_global_offset = [0.0]
    else:
        args_global_offset = [float(dum) for dum in
                              str(args.global_offset).split(",")]

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
    bigdict = json.loads(open(args.traces_file.name).read())

    # Load metadata from the traces
    meta_info = bigdict['meta_info']

    origin = meta_info['origin']
    insconf_uuid = origin['insconf_uuid']
    date_obs = origin.get('date_obs')

    tags = bigdict['tags']
    insmode = tags['insmode']

    # create instrument model
    pkg_paths = ['megaradrp.instrument.configs']
    store = asb.load_paths_store(pkg_paths)

    insmodel = asb.assembly_instrument(store, insconf_uuid, date_obs, by_key='uuid')

    pseudo_slit_config = insmodel.get_value('pseudoslit.boxes', **tags)

    fibid_with_box = assign_boxes_to_fibers(pseudo_slit_config, insmode)
    total_fibers = bigdict['total_fibers']
    if total_fibers != len(fibid_with_box):
        raise ValueError('Mismatch between number of fibers and '
                         'expected number from account from boxes')
    if 'global_offset' in bigdict.keys():
        global_offset = bigdict['global_offset']
        if args_global_offset != [0.0] and global_offset != [0.0]:
            raise ValueError('global_offset != 0 argument cannot be employed '
                             'when global_offset != 0 in JSON file')
        elif args_global_offset != [0.0]:
            global_offset = args_global_offset
    else:
        global_offset = args_global_offset
    print('>>> Using global_offset:', global_offset)
    pol_global_offset = np.polynomial.Polynomial(global_offset)
    if 'ref_column' in bigdict.keys():
        ref_column = bigdict['ref_column']
    else:
        ref_column = 2000
    for fiberdict in bigdict['contents']:
        fibid = fiberdict['fibid']
        fiblabel = fibid_with_box[fibid - 1]
        start = fiberdict['start']
        stop = fiberdict['stop']
        coeff = np.array(fiberdict['fitparms'])
        # skip fibers without trace
        if len(coeff) > 0:
            pol_trace = np.polynomial.Polynomial(coeff)
            y_at_ref_column = pol_trace(ref_column)
            correction = pol_global_offset(y_at_ref_column)
            coeff[0] += correction
            # update values in bigdict (JSON structure)
            bigdict['contents'][fibid-1]['fitparms'] = coeff.tolist()
            plot_trace(ax, coeff, start, stop, ix_offset, args.rawimage,
                       args.fibids, fiblabel, colour='blue')
        else:
            print('Warning ---> Missing fiber:', fibid_with_box[fibid - 1])

    # if present, read healing JSON file
    if args.healing is not None:
        healdict = json.loads(open(args.healing.name).read())
        list_operations = healdict['operations']
        for operation in list_operations:

            if operation['description'] == 'vertical_shift_in_pixels':
                if 'fibid_list' in operation.keys():
                    fibid_list = operation['fibid_list']
                else:
                    fibid_ini = operation['fibid_ini']
                    fibid_end = operation['fibid_end']
                    fibid_list = range(fibid_ini, fibid_end + 1)
                for fibid in fibid_list:
                    if fibid < 1 or fibid > total_fibers:
                        raise ValueError('fibid number outside valid range')
                    fiblabel = fibid_with_box[fibid - 1]
                    coeff = np.array(
                        bigdict['contents'][fibid - 1]['fitparms']
                    )
                    if len(coeff) > 0:
                        if args.verbose:
                            print('(vertical_shift_in_pixels) fibid:',
                                  fiblabel)
                        vshift = operation['vshift']
                        coeff[0] += vshift
                        bigdict['contents'][fibid - 1]['fitparms'] = \
                            coeff.tolist()
                        start = bigdict['contents'][fibid - 1]['start']
                        stop = bigdict['contents'][fibid - 1]['stop']
                        plot_trace(ax, coeff, start, stop, ix_offset,
                                   args.rawimage, True, fiblabel,
                                   colour='green')
                    else:
                        print('(vertical_shift_in_pixels SKIPPED) fibid:',
                              fiblabel)

            elif operation['description'] == 'duplicate_trace':
                fibid_original = operation['fibid_original']
                if fibid_original < 1 or fibid_original > total_fibers:
                    raise ValueError(
                        'fibid_original number outside valid range'
                    )
                fibid_duplicated = operation['fibid_duplicated']
                if fibid_duplicated < 1 or fibid_duplicated > total_fibers:
                    raise ValueError(
                        'fibid_duplicated number outside valid range'
                    )
                fiblabel_original = fibid_with_box[fibid_original - 1]
                fiblabel_duplicated = fibid_with_box[fibid_duplicated - 1]
                coeff = np.array(
                    bigdict['contents'][fibid_original - 1]['fitparms']
                )
                if len(coeff) > 0:
                    if args.verbose:
                        print('(duplicated_trace) fibids:',
                              fiblabel_original, '-->', fiblabel_duplicated)
                    vshift = operation['vshift']
                    coeff[0] += vshift
                    bigdict['contents'][fibid_duplicated - 1]['fitparms'] = \
                        coeff.tolist()
                    start = bigdict['contents'][fibid_original - 1]['start']
                    stop = bigdict['contents'][fibid_original - 1]['stop']
                    bigdict['contents'][fibid_duplicated - 1]['start'] = start
                    bigdict['contents'][fibid_duplicated - 1]['stop'] = stop
                    plot_trace(ax, coeff, start, stop, ix_offset,
                               args.rawimage, True, fiblabel_duplicated,
                               colour='green')
                else:
                    print('(duplicated_trace SKIPPED) fibids:',
                          fiblabel_original, '-->', fiblabel_duplicated)

            elif operation['description'] == 'extrapolation':
                if 'fibid_list' in operation.keys():
                    fibid_list = operation['fibid_list']
                else:
                    fibid_ini = operation['fibid_ini']
                    fibid_end = operation['fibid_end']
                    fibid_list = range(fibid_ini, fibid_end + 1)
                for fibid in fibid_list:
                    if fibid < 1 or fibid > total_fibers:
                        raise ValueError('fibid number outside valid range')
                    fiblabel = fibid_with_box[fibid - 1]
                    coeff = np.array(
                        bigdict['contents'][fibid - 1]['fitparms']
                    )
                    if len(coeff) > 0:
                        if args.verbose:
                            print('(extrapolation) fibid:', fiblabel)
                        # update values in bigdict (JSON structure)
                        start = operation['start']
                        stop = operation['stop']
                        start_orig = bigdict['contents'][fibid - 1]['start']
                        stop_orig = bigdict['contents'][fibid - 1]['stop']
                        bigdict['contents'][fibid - 1]['start'] = start
                        bigdict['contents'][fibid - 1]['stop'] = stop
                        if start < start_orig:
                            plot_trace(ax, coeff, start, start_orig,
                                       ix_offset,
                                       args.rawimage, True, fiblabel,
                                       colour='green')
                        if stop_orig < stop:
                            plot_trace(ax, coeff, stop_orig, stop,
                                       ix_offset,
                                       args.rawimage, True, fiblabel,
                                       colour='green')
                        if start_orig <= start <= stop <= stop_orig:
                            plot_trace(ax, coeff, start, stop,
                                       ix_offset,
                                       args.rawimage, True, fiblabel,
                                       colour='green')
                    else:
                        print('(extrapolation SKIPPED) fibid:', fiblabel)

            elif operation['description'] == 'fit_through_user_points':
                fibid = operation['fibid']
                fiblabel = fibid_with_box[fibid - 1]
                if args.verbose:
                    print('(fit through user points) fibid:', fiblabel)
                poldeg = operation['poldeg']
                start = operation['start']
                stop = operation['stop']
                xfit = []
                yfit = []
                for userpoint in operation['user_points']:
                    # assume x, y coordinates in JSON file are given in
                    # image coordinates, starting at (1,1) in the lower
                    # left corner
                    xdum = userpoint['x'] - 1  # use np.array coordinates
                    ydum = userpoint['y'] - 1  # use np.array coordinates
                    xfit.append(xdum)
                    yfit.append(ydum)
                xfit = np.array(xfit)
                yfit = np.array(yfit)
                if len(xfit) <= poldeg:
                    raise ValueError('Insufficient number of points to fit'
                                     ' polynomial')
                poly, residum = polfit_residuals(xfit, yfit, poldeg)
                coeff = poly.coef
                plot_trace(ax, coeff, start, stop, ix_offset,
                           args.rawimage, args.fibids, fiblabel,
                           colour='green')
                bigdict['contents'][fibid - 1]['start'] = start
                bigdict['contents'][fibid - 1]['stop'] = stop
                bigdict['contents'][fibid - 1]['fitparms'] = coeff.tolist()

            elif operation['description'] == \
                    'extrapolation_through_user_points':
                fibid = operation['fibid']
                fiblabel = fibid_with_box[fibid - 1]
                if args.verbose:
                    print('(extrapolation_through_user_points):', fiblabel)
                start_reuse = operation['start_reuse']
                stop_reuse = operation['stop_reuse']
                resampling = operation['resampling']
                poldeg = operation['poldeg']
                start = operation['start']
                stop = operation['stop']
                coeff = bigdict['contents'][fibid - 1]['fitparms']
                xfit = np.linspace(start_reuse, stop_reuse, num=resampling)
                poly = np.polynomial.Polynomial(coeff)
                yfit = poly(xfit)
                for userpoint in operation['user_points']:
                    # assume x, y coordinates in JSON file are given in
                    # image coordinates, starting at (1,1) in the lower
                    # left corner
                    xdum = userpoint['x'] - 1  # use np.array coordinates
                    ydum = userpoint['y'] - 1  # use np.array coordinates
                    xfit = np.concatenate((xfit, np.array([xdum])))
                    yfit = np.concatenate((yfit, np.array([ydum])))
                poly, residum = polfit_residuals(xfit, yfit, poldeg)
                coeff = poly.coef
                if start < start_reuse:
                    plot_trace(ax, coeff, start, start_reuse, ix_offset,
                               args.rawimage, args.fibids, fiblabel,
                               colour='green')
                if stop_reuse < stop:
                    plot_trace(ax, coeff, stop_reuse, stop, ix_offset,
                               args.rawimage, args.fibids, fiblabel,
                               colour='green')
                bigdict['contents'][fibid - 1]['start'] = start
                bigdict['contents'][fibid - 1]['stop'] = stop
                bigdict['contents'][fibid - 1]['fitparms'] = coeff.tolist()

            elif operation['description'] == 'sandwich':
                fibid = operation['fibid']
                fiblabel = fibid_with_box[fibid - 1]
                if args.verbose:
                    print('(sandwich) fibid:', fiblabel)
                fraction = operation['fraction']
                nf1, nf2 = operation['neighbours']
                start = operation['start']
                stop = operation['stop']
                tmpf1 = bigdict['contents'][nf1 - 1]
                tmpf2 = bigdict['contents'][nf2 - 1]
                if nf1 != tmpf1['fibid'] or nf2 != tmpf2['fibid']:
                    raise ValueError(
                        "Unexpected fiber numbers in neighbours"
                    )
                coefff1 = np.array(tmpf1['fitparms'])
                coefff2 = np.array(tmpf2['fitparms'])
                coeff = coefff1 + fraction * (coefff2 - coefff1)
                plot_trace(ax, coeff, start, stop, ix_offset,
                           args.rawimage, args.fibids,
                           fiblabel, colour='green')
                # update values in bigdict (JSON structure)
                bigdict['contents'][fibid - 1]['start'] = start
                bigdict['contents'][fibid - 1]['stop'] = stop
                bigdict['contents'][fibid - 1][
                    'fitparms'] = coeff.tolist()
                if fibid in bigdict['error_fitting']:
                    bigdict['error_fitting'].remove(fibid)

            elif operation['description'] == 'renumber_fibids_within_box':
                fibid_ini = operation['fibid_ini']
                fibid_end = operation['fibid_end']
                box_ini = fibid_with_box[fibid_ini - 1][4:]
                box_end = fibid_with_box[fibid_end - 1][4:]
                if box_ini != box_end:
                    print(f'ERROR: box_ini={box_ini}, box_end={box_end}')
                    raise ValueError('fibid_ini and fibid_end correspond to '
                                     'different fiber boxes')
                fibid_shift = operation['fibid_shift']
                if fibid_shift in [-1, 1]:
                    if fibid_shift == -1:
                        i_start = fibid_ini
                        i_stop = fibid_end + 1
                        i_step = 1
                    else:
                        i_start = fibid_end
                        i_stop = fibid_ini - 1
                        i_step = -1
                    for fibid in range(i_start, i_stop, i_step):
                        fiblabel_ori = fibid_with_box[fibid - 1]
                        fiblabel_new = fibid_with_box[fibid - 1 + fibid_shift]
                        if args.verbose:
                            print('(renumber_fibids) fibid:',
                                  fiblabel_ori, '-->', fiblabel_new)
                        bigdict['contents'][fibid -1 + fibid_shift] = \
                            deepcopy(bigdict['contents'][fibid -1])
                        bigdict['contents'][fibid -1 + fibid_shift]['fibid'] += \
                            fibid_shift
                        # display updated trace
                        coeff = \
                            bigdict['contents'][fibid -1 + fibid_shift]['fitparms']
                        start = \
                            bigdict['contents'][fibid -1 + fibid_shift]['start']
                        stop = bigdict['contents'][fibid -1 + fibid_shift]['stop']
                        plot_trace(ax, coeff, start, stop, ix_offset,
                                   args.rawimage, args.fibids,
                                   fiblabel_ori + '-->' + fiblabel_new,
                                   colour='green')
                    if fibid_shift == -1:
                        bigdict['contents'][fibid_end - 1]['fitparms'] = []
                    else:
                        bigdict['contents'][fibid_ini - 1]['fitparms'] = []
                else:
                    raise ValueError('fibid_shift in operation '
                                     'renumber_fibids_within_box '
                                     'must be -1 or 1')
            else:
                raise ValueError('Unexpected healing method:',
                                 operation['description'])

# update trace map
    if args.updated_traces is not None:
        # avoid overwritting initial JSON file
        if args.updated_traces.name != args.traces_file.name:
            # new random uuid for the updated calibration
            bigdict['uuid'] = str(uuid4())
            with open(args.updated_traces.name, 'w') as outfile:
                json.dump(bigdict, outfile, indent=2)

    if pdf is not None:
        pdf.savefig()
        pdf.close()
    else:
        pause_debugplot(12, pltshow=True, tight_layout=True)


if __name__ == "__main__":

    main()
