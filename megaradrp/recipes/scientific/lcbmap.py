import math
import random
import logging
import numpy as np
import matplotlib.pylab as plt
from scipy.spatial import KDTree

from abc import ABCMeta, abstractmethod

logger = logging.getLogger(__name__)


class Map(object):
    """
    An top level object for managing all game data related to positioning.
    """
    directions = [(0, 1), (1, 1), (1, 0), (0, -1), (-1, -1), (-1, 0)]

    def __init__(self, ( rows, cols ), *args, **keywords):
        # Map size
        self.rows = rows
        self.cols = cols

    def __str__(self):
        return "Map (%d, %d)" % (self.rows, self.cols)

    @property
    def size(self):
        """Returns the size of the grid as a tuple (row, col)"""
        return (self.rows, self.cols)

    def distance(self, start, destination):
        """Takes two hex coordinates and determine the distance between them."""
        logger.debug("Start: %s, Dest: %s", start, destination)
        diffX = destination[0] - start[0]
        diffY = destination[1] - start[1]

        distance = min(abs(diffX), abs(diffY)) + abs(diffX - diffY)

        logger.debug("diffX: %d, diffY: %d, distance: %d", diffX, diffY,
                     distance)
        return distance

    @classmethod
    def direction(self, origin, destination):
        """
        Reports the dominating direction from an origin to a destination.  if even, chooses randomly
        Useful for calculating any type of forced movement
        """
        offset = (destination[0] - origin[0], destination[1] - origin[1])
        scale = float(max(abs(offset[0]), abs(offset[1])))
        if scale == 0:
            return (0, 0)
        direction = (offset[0] / scale, offset[1] / scale)

        # Handle special cases
        if direction == (1, -1):
            return random.choice([(1, 0), (0, -1)])
        elif direction == (-1, 1):
            return random.choice([(-1, 0), (0, 1)])

        def choose(i):
            if i == 0.5:
                return random.choice((0, 1))
            elif i == -0.5:
                return random.choice((0, -1))
            else:
                return int(round(i))

        return (choose(direction[0]), choose(direction[1]))

    def ascii(self, numbers=True):
        """ Debug method that draws the grid using ascii text """

        table = ""

        if numbers:
            text_length = len(str(self.rows - 1 if self.cols % 2 == 1 else self.rows) + ',' + str(int(self.rows - 1 + math.floor(self.cols / 2))))
        else:
            text_length = 3

            # Header for first row
        for col in range(self.cols):
            if col % 2 == 0:
                table += " " + '_' * text_length
            else:
                table += " " + ' ' * text_length
        table += "\n"
        # Each row
        for row in range(self.rows):
            top = "/"
            bottom = "\\"

            for col in range(self.cols):
                if col % 2 == 0:
                    # text = "%d,%d" % (row + col / 2, col) if numbers else ""
                    text = "%d" % (self.units.default[row, col]) if numbers and self.units.default[row, col] else ""
                    top += (text).center(text_length) + "\\"
                    bottom += ("").center(text_length, '_') + "/"
                else:
                    text = "%d" % (self.units.default[row, col]) if numbers and self.units.default[row, col] else " "
                    top += ("").center(text_length, '_') + "/"
                    bottom += (text).center(text_length) + "\\"
            # Clean up tail slashes on even numbers of columns
            if self.cols % 2 == 0:
                if row == 0:
                    top = top[:-1]
            table += top + "\n" + bottom + "\n"

        # Footer for last row
        footer = " "
        for col in range(0, self.cols - 1, 2):
            footer += " " * text_length + "\\" + '_' * text_length + "/"
        table += footer + "\n"
        return table

    def valid_cell(self, cell):
        row, col = cell
        if col < 0 or col >= self.cols: return False
        if row < math.ceil(col / 2.0) or row >= math.ceil(
                        col / 2.0) + self.rows: return False
        return True

    def neighbors(self, center):
        """
        Return the valid cells neighboring the provided cell.
        """
        return filter(self.valid_cell, [
            (center[0] - 1, center[1]), (center[0], center[1] + 1),
            (center[0] + 1, center[1] + 1), (center[0] + 1, center[1]),
            (center[0], center[1] - 1), (center[0] - 1, center[1] - 1)
        ])

    def spread(self, center, radius=1):
        """
        A slice of a map is a collection of valid cells, starting at an origin,
        and encompassing all cells within a given radius.
        """
        result = set((center,))  # Start out with this center cell
        neighbors = self.neighbors(center)  # Get the neighbors for use later
        if radius == 1:  # Recursion end case
            result = result | set(
                neighbors)  # Return the set of this cell and its neighbors
        else:  # Otherwise, recurse over all the neghbors,
            for n in neighbors:  # decrementing the radius by one.
                result = result | set(self.spread(n, radius - 1))
        return filter(self.valid_cell,
                      result)  # filter invalid cells before returning.

    def cone(self, origin, direction, length=1):
        """
        A slice of a map is a section of cells, originating a a single cell and
        extending outward through two cells separated by one cell between them.
        In the example below, starting at (0,0), (0,1) and (1,0) define a slice,
        as do (-1,-1) and (0,1).
               _____
         _____/-1,0 \_____
        /-1,-1\_____/ 0,1 \
        \_____/ 0,0 \_____/
        /0,-1 \_____/ 1,1 \
        \_____/ 1,0 \_____/
              \_____/
        """
        result = self.slice(origin, direction, length)
        result.extend(self.slice(origin, (direction + 1) % 6, length))
        return filter(self.valid_cell, set(result))

    def slice(self, origin, direction, length=2):
        """
        A slice of a map is a section of cells, originating a a single cell and
        extending outward through two neighboring cells.  In the example below,
        starting at (0,0), (0,1) and (1,1) define a slice, as do (-1,0) and
        (-1,-1).
               _____
         _____/-1,0 \_____
        /-1,-1\_____/ 0,1 \
        \_____/ 0,0 \_____/
        /0,-1 \_____/ 1,1 \
        \_____/ 1,0 \_____/
              \_____/
        """
        # The edge wheel described in the docnotes above, used for calculating edges and steps


        # edge is the step we take for each distance,
        # step is the increment for each cell that distance out
        edge, step = self.directions[direction % 6], self.directions[
            (direction + 2) % 6]

        logger.debug("Edge: %s, Step: %s", edge, step)

        result = [origin]
        # Work each row, i units out along an edge
        for i in range(1, length + 1):
            start = (origin[0] + edge[0] * i, origin[1] + edge[1] * i)
            for j in range(i + 1):
                # calculate
                pos = (start[0] + step[0] * j, start[1] + step[1] * j)
                result.append(pos)
        return filter(self.valid_cell, result)

    def line(self, origin, direction, length=3):
        """
        Returns all the cells along a given line, starting at an origin
        """
        offset = self.directions[direction]
        results = [origin]
        # Work each row, i units out along an edge
        for i in range(1, length + 1):
            results.append(
                (origin[0] + offset[0] * i, origin[1] + offset[1] * i))
        return filter(self.valid_cell, results)

    def cells(self):
        cells = []
        for row in range(self.rows + self.cols / 2):
            cells.extend(
                [(row, col) for col in range(1 + 2 * row)]
            )
        return filter(self.valid_cell, cells)


class Grid(dict):
    """An extension of a basic dictionary with a fast, consistent lookup by value implementation."""

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
        print df.to_string()

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
                fibra = value['real_fiber']
                x_pos = value['x_pos'][indice]
                y_pos = value['y_pos'][indice]
                return pos, indice, trace, fibra, x_pos, y_pos
        return None

    def get_neighbours(self, points, radius=0.5):
        ind = np.array(self.neighbours.query_ball_point(points, r=radius))+1 # We add +1 because .dat file starts in 1
        return ind

class MapUnit(object):
    """
    An abstract base class that will contain or require implementation of all
    the methods necessary for a unit to be managed by a map object.
    """

    __metaclass__ = ABCMeta

    def __init__(self, grid):
        self.grid = grid

    @property
    def position(self):
        """A property that looks up the position of this unit on it's associated map."""
        return self.grid.find(self)

    @abstractmethod
    def paint(self, surface):
        """An abstract base method to contain the painting code for a given unit."""
        pass



if __name__ == '__main__':

    m = Map((21, 27))
    m.units = Grid('LCB_spaxel_centers_new.dat')
    numbers = True

    test_cell = (1, 1)

    print (m.ascii(numbers))
    # print ("Valid Cell?: %s" % m.valid_cell(test_cell))
    # print ("Neighbours: %s" % m.neighbors(test_cell))
    # print ("Object: %s" % m)
    # print ("Spread: %s" % m.spread(test_cell))
    # print ("line: %s" % m.line(test_cell, 0))

    # print ("distancia de (1,1) a (2,3): %s" %m.distance((1, 1), (2, 3)))

    print ("lista: ", [289, 290, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 334, 335])
    for trace in [289, 290, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 334, 335]:
        print ("Datos de la fibra %s (hueco, indice_diccionario, fibra_real, linea_en_fits, x_pos, y_pos): %s" %(trace,m.units.get_from_trace(trace)))

    points = [0,0]
    vecinos = m.units.get_neighbours(points, radius=0.5)
    for elem in vecinos:
        print 'vecinos', m.units.get_from_trace(elem)


    # spaxels = np.loadtxt('LCB_spaxel_centers_new.dat')
    # lines = spaxels[:,4].astype(int)
    # bundles = spaxels[:,5].astype(int)
    # fibers = spaxels[:,3].astype(int)
    # x_pos = spaxels[:,0]
    # y_pos = spaxels[:,1]

    # from scipy.spatial import KDTree
    #
    # data = zip(x_pos.ravel(), y_pos.ravel())
    # arbol = KDTree(data)
    # points = [[0,0]]
    # distances, indices = arbol.query(points, k=5)
    # print distances
    # print indices
    # print(repr(indices))
    #
    # ind = arbol.query_ball_point(points, r=0.5)
    # print ind