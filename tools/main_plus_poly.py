from __future__ import print_function

import glob
import shutil
import json

import pandas as pd

import os


def copy_from_folder():
    for base in glob.glob("*"):
        print(base)
        for file in glob.glob("%s/*.json" % base):
            full_file_name = "%s" % file
            dest = "%s.json" % base
            print("\t %s" % full_file_name)

            shutil.copy(full_file_name, dest)


def read_json(name):
    lin_lrv = [6163.59390, 6150.29600, 6143.06260, 6128.44990, 6096.16310,
               6074.33770, 6029.99690, 5975.53400, 5944.83420, 5913.63100,
               5881.89520, 5852.48790, 5433.64990, 5418.55770, 5400.56180,
               5341.09320]

    lin_lri = [7272.935, 7372.128, 7383.980, 7435.488, 7503.868, 7514.651,
               7635.105, 7723.984, 7891.040, 7948.176, 7979.004, 8006.156,
               8014.785, 8053.307, 8103.692, 8115.311, 8264.521, 8330.425,
               8408.209, 8424.647, 8521.441, 8605.768, 8620.491, 8667.943]

    lin_lrz = [8046.11, 8136.40, 8205.11, 8266.08, 8320.86, 8358.72, 8421.22,
               8446.51, 8582.90, 8621.30, 8709.23, 8812.51, 8865.31, 9148.67,
               9276.27, 9459.21, 9547.74]

    lin_lrr = [7298.1436, 7173.3726, 7059.5254, 6955.3149, 6876.2925,
               6658.6772, 6524.7627, 6412.5918, 6337.6206, 6326.3667,
               6214.4413, 6098.8032, 6090.1079, 6020.1191, 6013.2817]

    with open('%s.json' % name) as data_file:
        data = json.load(data_file)

    c0 = []
    c1 = []
    c2 = []
    c3 = []
    c4 = []
    c5 = []

    xpos = []
    ypos = []
    fwhm = []
    wave = []
    ref = []
    list_fib = []
    list_coeffs = []
    list_ind = []

    fibra = 1

    contents = data['contents']
    if isinstance(contents, dict):
        conts = contents.values()
    else:
        conts = contents
    for elem in conts:
        list_coeffs=elem['coeff']

        for linea in elem['features']:
            xpos.append(linea['xpos'])
            ypos.append(linea['ypos'])
            fwhm.append(linea['fwhm'])
            wave.append(linea['wavelength'])
            ref.append(linea['reference'])
            c0.append(list_coeffs[0])
            c1.append(list_coeffs[1])
            c2.append(list_coeffs[2])
            c3.append(list_coeffs[3])
            c4.append(list_coeffs[4])
            c5.append(list_coeffs[5])
            list_fib.append(fibra)
            if 'LRZ' in name:
                list_ind.append(lin_lrz.index(linea['reference']) + 1)
            elif 'LRV' in name:
                list_ind.append(lin_lrv.index(linea['reference']) + 1)
            elif 'LRR' in name:
                try:
                    list_ind.append(lin_lrr.index(linea['reference']) + 1)
                except:
                    list_ind.append(0)
            else:
                try:
                    list_ind.append(lin_lri.index(linea['reference']) + 1)
                except:
                    list_ind.append(0)
        fibra += 1

    data = {
            'xpos': xpos,
            'ypos': ypos,
            'fwhm': fwhm,
            'wavelength': wave,
            'reference': ref,
            'fibra': list_fib,
            'line': list_ind,
            'zc0': c0,
            'zc1': c1,
            'zc2': c2,
            'zc3': c3,
            'zc4': c4,
            'zc5': c5
            }
    df = pd.DataFrame(data)
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    excel = pd.ExcelWriter('%s.xlsx' % name, engine='xlsxwriter')
    df.to_excel(excel, sheet_name='Sheet1', index=False)
    # df.to_excel(writer, sheet_name='Sheet2')
    csv = df.to_csv('%s.csv' % name, sep=' ', index=False, encoding='utf-8')


def generar_xls():
    lista_ficheros = []
    for base in glob.glob("*.json"):
        file_name = base.split('.')[0]
        lista_ficheros.append(file_name)
        read_json(file_name)


    lista = {'LRR':[],
             'LRV':[],
             'LRZ':[],
             'sci_LRR':[],
             'sci_LRV':[],
             'sci_LRZ':[],
             'sci_LRI':[],
             }
    for elem in lista_ficheros:
        if 'LRZ' in elem:
            if 'sci' in elem:
                lista['sci_LRZ'].append(elem)
            else:
                lista['LRZ'].append(elem)
        elif 'LRR' in elem:
            if 'sci' in elem:
                lista['sci_LRR'].append(elem)
            else:
                lista['LRR'].append(elem)
        elif 'LRI' in elem:
            if 'sci' in elem:
                lista['sci_LRI'].append(elem)
            else:
                lista['LRI'].append(elem)
        else:
            if 'sci' in elem:
                lista['sci_LRV'].append(elem)
            else:
                lista['LRV'].append(elem)

    for elem in lista:
        lista[elem].sort()
        if lista[elem]:
            file = open('%s.txt'%elem, 'w')
            for item in lista[elem]:
                file.write("%s\n" % item)


os.chdir(".")
generar_xls()

print('**************************Fin*************************')
