
from __future__ import print_function

import logging
import numpy as np
from scipy.spatial import KDTree

logger = logging.getLogger(__name__)

class Grid(dict):
    """
    """

    def __init__(self, spaxels, *args, **keywords):

        lines = spaxels[:,3].astype(int)
        fibers = spaxels[:,2].astype(int)
        x_pos = spaxels[:,0]
        y_pos = spaxels[:,1]

        neighbours = zip(x_pos.ravel(), y_pos.ravel())
        self.neighbours = KDTree(neighbours)

        diccionario= dict((k, []) for k in np.unique(fibers))

        for fib in np.unique(fibers):
            aux = {}
            aux['real_fiber'] = fib
            aux['trace'] = lines[fibers==fib].tolist()
            aux['x_pos'] = x_pos[fibers==fib].tolist()
            aux['y_pos'] = y_pos[fibers==fib].tolist()

            diccionario[fib] = aux

        super(Grid, self).__init__(diccionario)

        # matriz = self.matrix_generation(spaxels)
        # matriz = self.delete_ceroes(matriz)
        #
        # matriz = np.flipud(matriz)
        #
        # self.default = matriz


    def __getitem__(self, key):
        return super(Grid, self).get(key, self.default)

    def matrix_generation(self, spaxels):
        spaxels = spaxels[np.where((spaxels[:,0] >= -20) & (spaxels[:,0] <=20) & (spaxels[:,1] >= -20) & (spaxels[:,1] <=20))]
        ox = np.around(spaxels[:,0], decimals=3)
        oy = np.around(spaxels[:,1], decimals=3)
        fiber_number = spaxels[:,3].astype(int)

        paso_ox = 0.465
        paso_oy = 0.268

        col_x = np.around(ox/paso_ox, decimals=1) + 13
        col_y = np.around(oy/paso_oy, decimals=1) + 21  #20 si usamos LCB_spaxel_centers.dat

        col_x = col_x.astype(int)
        col_y = col_y.astype(int)

        coordenadas = zip(col_x, col_y, fiber_number)

        matriz = np.zeros((42,27))

        for elem in coordenadas:
            matriz[elem[1],elem[0]] = elem[2]

        return matriz

    def delete_ceroes(self, matriz):
        # pintar_matriz(matriz)
        nueva_matriz = np.zeros((22,27))
        for indice, row in enumerate(matriz):
            nueva_fila = row
            if indice not in [matriz.shape[0]-1,0]:
                if indice%2 ==0:
                    nueva_fila = row + matriz[indice-1,:]
                    nueva_matriz[indice/2,:] = nueva_fila
                else:
                    nueva_matriz[(indice//2)+1,:] = nueva_fila
            elif indice in [0]:
                nueva_matriz[indice,:] = nueva_fila
            else:
                nueva_matriz[-1,:] = nueva_fila
        for cont in xrange(1,27,2):
            nueva_matriz[:,cont] = np.roll(nueva_matriz[:,cont].T, 1, axis=0).T
        # nueva_matriz = np.delete(nueva_matriz, -1, 0)

        return nueva_matriz

    def dump_matriz(self, matriz):
        import pandas
        df = pandas.DataFrame(matriz) # Papel
        # df = pandas.DataFrame(matriz)
        print(df.to_string())

    def find(self, item):
        """
        A fast lookup by value implementation
        """
        for pos, value in self.items():
            if item == value or item in value:
                return pos
        return None

    def get_from_trace(self, trace):
        for pos, value in self.items():
            if trace in value['trace']:
                indice = np.where(np.in1d(value['trace'], trace)==True)[0][0]
                bun = value['real_fiber']
                x_pos = value['x_pos'][indice]
                y_pos = value['y_pos'][indice]
                return trace, bun, x_pos, y_pos
        return None

    def get_neighbours(self, points, radius=0.5):
        ind = np.array(self.neighbours.query_ball_point(points, r=radius))+1 # We add +1 because .dat file starts in 1
        return ind

    def get_fiber(self, points):
        ind = np.array(self.neighbours.query(points, k=1))[1]
        return int(ind)