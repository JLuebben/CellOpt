from __future__ import print_function
from copy import deepcopy
from math import cos, pi
from subprocess import call, STDOUT
from shutil import copyfile
import os
import sys
import argparse
from os.path import dirname, join
from collections import OrderedDict
try:
    import matplotlib.pyplot as plt
    PLOTAVAILABLE = True
except ImportError:
    PLOTAVAILABLE = False
    plt = None


CLASSPARAMETERS = {'triclinic': ((0, 1, 2, 3, 4, 5), {}),
                   'monoclinic': ((0, 1, 2, 4), {}),
                   'orthorhombic': ((0, 1, 2), {}),
                   'tetragonal': ((0, 2), {0: (1,)}),
                   'rhombohedral': ((0, 3), {0: (1, 2), 3: (4, 5)}),
                   'hexagonal': ((0, 2), {0: (1,)}),
                   'cubic': ((0,), {0: (1, 2)})}


def callShelxl(fileName):
    """
    Call SHELXL in a subprocess.
    :param fileName: str
    :return: None
    """
    FNULL = open(os.devnull, 'w')
    try:
        call(['shelxl.exe', fileName], stdout=FNULL, stderr=STDOUT)
    except:
        call(['shelxl', fileName], stdout=FNULL, stderr=STDOUT)


def evaluate(fileName):
    """
    Call SHELXL and subsequently evaluate the result.
    :param fileName: str
    :return: float<wR2>, float<meanDfixFit>, float<weightedDfixFit>
    """
    callShelxl(fileName)
    wR2 = 999
    with open('work.lst', 'r') as fp:
        for line in fp.readlines():
            if 'for all data' in line:
                line = [word for word in line.split() if line]
                wR2 = float(line[2][:-1])
                break
    reader = ShelxlReader()
    molecule = reader.read(fileName + '.res')
    mean, weighted = None, None
    try:
        mean, weighted = molecule.checkDfix()
    except ZeroDivisionError:
        print('\n\n\nSomething went wrong while re-refining the structure.')
        print('\n\nError Messages from work.lst file:')
        with open('work.lst', 'r') as fp:
            for line in fp.readlines():
                if '**' in line:
                    print(line[:-1])
        print('\nExiting')
        exit(1)
    except ValueError:
        print('\n\n\nSomething went wrong while re-refining the structure.')
        print('\n\nError Messages from work.lst file:')
        with open('work.lst', 'r') as fp:
            for line in fp.readlines():
                if '**' in line:
                    print(line[:-1])
                    break
        print('\nExiting')
        exit(1)
        # print('\n\n\nNo DFIX or DANG restraints found in structure.\n\nExiting')
        # exit(2)
    return wR2, mean, weighted


def quickEvaluate(molecule, cell):
    """
    Evaluates the DFIX fit to a given cell
    :param molecule: ShelxlMolecule instance
    :param cell: list of six floats
    :return: float<meanDfixFit>, float<weightedDfixFit>
    """
    molecule.cell = cell
    try:
        return molecule.checkDfix()
    except ValueError:
        print('\n\n\nNo DFIX or DANG restraints found in structure.\n\nExiting')
        exit(2)


def determineCrystalClass(cell):
    """
    Derive crystal class from the cell parameter values of a given cell
    :param cell: list of six floats
    :return: str<className>, tuple<constraints>
    """
    a, b, c, A, B, C = cell[2:]
    cls = 'triclinic'
    if a == b == c and A == B == C and int(float(A)) == 90:
        cls = 'cubic'
    elif a == b and A == B == C and int(float(A)) == 120:
        cls = 'hexagonal'
    elif a == b == c and A == B == C and not int(float(A)) == 90:
        cls = 'rhombohedral'
    elif a == b and A == B == C and int(float(A)) == 90:
        cls = 'tetragonal'
    elif a != b != c and A == B == C and int(float(A)) == 90:
        cls = 'orthorhombic'
    elif a != b != c and A == C and int(float(A)) == 90 and not int(float(B)) == 90:
        cls = 'monoclinic'
    return cls, CLASSPARAMETERS[cls]


def generateJobs(params, cell, delta):
    """
    Generate optimization job steps based on a tuple defining constraints, unit cell parameters and the current step
    size
    :param params: tuple<constraints>
    :param cell: list of six floats
    :param delta: float
    :return: list of new cells
    """
    data = [float(x) for x in cell[2:]]
    jobs = [data]
    refine, conDict = params
    for p in params[0]:
        j = data[:]
        j[p] -= delta
        for source, targets in conDict.items():
            for target in targets:
                j[target] = j[source]
        jobs.append(j)

        j = data[:]
        j[p] += delta
        for source, targets in conDict.items():
            for target in targets:
                j[target] = j[source]
        jobs.append(j)
    return jobs


def run(fileName, p1=False, overrideClass=None, fast=False, plot=False):
    """
    Run the optimizer in 'fast' or 'default' mode.
    :param fileName: str<Name of the starting parameter shelxl.res file>
    :param p1: bool<Expand structure to P1/P-1>
    :param overrideClass: str<name of crystal class>
    :param fast: bool<use fast optimization scheme>
    :param plot: bool<plot diagnostics plot.>
    :return: None
    """
    plotter = Plotter()
    resFileName = fileName + '.res'
    fileDir = dirname(resFileName)
    copyfile(join(fileDir, fileName + '.hkl'), './work.hkl')
    reader = ShelxlReader()
    molecule = reader.read(resFileName)
    cell = reader['cell'].split()

    cls, params = determineCrystalClass(cell)
    if overrideClass:
        cls = overrideClass
        params = CLASSPARAMETERS[cls]
    print('Crystal Class is {}.'.format(cls))
    if p1:
        if '-' in reader['latt']:
            print('Expanding to P1.')
        else:
            print('Expanding to P-1.')
        try:
            reader.toP1()
        except ValueError:
            print('Expanding Structures to P1/P-1 is not supported for structures\n'
                  'containing multiple residues.')
            exit(3)
        cls = 'triclinic'
        params = CLASSPARAMETERS[cls]
    originalCell = [float(x) for x in cell[2:]]
    startDiff = 999
    try:
        startDiff, _ = molecule.checkDfix()
    except ValueError:
        print('\nNo DFIX or DANG restraints found in structure.\n\nExiting')
        exit(2)
    startDiff0 = startDiff
    lastDiff = 9999

    iterations = 25

    i = -1
    barLengths = 20

    print('  ' + (barLengths - 8) // 2 * '-' + 'Progress' + (
                barLengths - 8) // 2 * '-' + '  ---Fit--   ---a---   ---b---   ---c---   -alpha-   --beta-   -gamma-')
    progress = (i + 1) / iterations
    progress = int(barLengths * progress)
    sys.stdout.write(
        '\r [' + progress * '#' + (barLengths - progress) * '-' + ']')
    sys.stdout.flush()
    for i in range(iterations):
        plotter(a=float(cell[2]), b=float(cell[3]), c=float(cell[4]), alpha=float(cell[5]), beta=float(cell[6]),
                gamma=float(cell[7]), fit=startDiff*100)
        i += 1
        sdelta = .1
        slastImprovement = 0
        for ii in range(250):
            sbestW = lastDiff
            sbestWj = 0
            jobs = generateJobs(params, cell, sdelta)
            for j, job in enumerate(jobs):
                weighted, mean = quickEvaluate(molecule, job)
                if weighted < sbestW:
                    sbestW = weighted
                    sbestWj = j
                    plotter(a=float(job[0]), b=float(job[1]), c=float(job[2]), alpha=float(job[3]), beta=float(job[4]),
                        gamma=float(job[5]), fit=sbestW*100)
                    progress = (i) / iterations
                    progress = int(barLengths * progress)
                    sys.stdout.write(
                        '\r [' + progress * '#' + (barLengths - progress) * '-' + '] {fit:8.6f} {cell}'.format(
                            fit=weighted,
                            cell=' '.join([
                                '{:9.4f}'.format(
                                    p)
                                for
                                p
                                in
                                job])))
                    sys.stdout.flush()
            cell = cell[:2] + ['{:7.4f}'.format(p) for p in jobs[sbestWj]]
            if sbestWj == 0:
                # print('No improvements found. Decreasing step size.')
                sdelta = sdelta / 2
                if sdelta < 0.002:
                    # print('Converged.')
                    break
                if ii - slastImprovement > 10:
                    # print('No improvements since 10 steps. Terminating.')
                    break
            else:
                slastImprovement = ii
        if not fast:
            newCell = ' '.join(cell) + '\n'
            reader['cell'] = newCell
            reader.write(fileName='work.ins')
            wR2, mean, weighted = evaluate('work')
            newReader = ShelxlReader()
            molecule = newReader.read('work.res')
            progress = i / iterations
            progress = int(barLengths * progress)
            sys.stdout.write(
                '\r [' + progress * '#' + (barLengths - progress) * '-'
                + '] {fit:8.6f} {cell}'.format(fit=weighted,
                                               cell=' '.join(['{:9.4f}'.format(p) for p in job])))
            sys.stdout.flush()
        else:
            progress = barLengths
            sys.stdout.write(
                '\r [' + progress * '#' + (barLengths - progress) * '-'
                + '] {fit:8.6f} {cell}'.format(fit=weighted,
                                               cell=' '.join(['{:9.4f}'.format(p) for p in job])))
            sys.stdout.flush()
            break
        startDiff = sbestW
    print()
    print('\n\nOriginal Cell:', cell2String(originalCell, offset=15))
    print()
    print('   Final Cell:', cell2String(cell[2:], offset=15))

    print('\nOriginal DFIX fit: {:8.6f}'.format(startDiff0))
    print('   Final DFIX fit: {:8.6f}'.format(weighted))
    if plot:
        plotter.show()


def run2(fileName, p1=False, overrideClass=None):
    """
    Run the optimizer in 'accurate' mode.
    :param fileName: str<Name of the starting parameter shelxl.res file>
    :param p1: bool<Expand structure to P1/P-1>
    :param overrideClass: str<name of crystal class>
    :return: None
    """
    resFileName = fileName + '.res'
    fileDir = dirname(resFileName)
    copyfile(join(fileDir, fileName + '.hkl'), './work.hkl')
    reader = ShelxlReader()
    molecule = reader.read(resFileName)
    cell = reader['cell'].split()
    cls, params = determineCrystalClass(cell)
    if overrideClass:
        cls = overrideClass
        params = CLASSPARAMETERS[cls]
    print('Crystal Class is {}.'.format(cls))
    if p1:
        print('Expanding to P1.')
        reader.toP1()
        cls = 'triclinic'
        params = CLASSPARAMETERS[cls]
    originalCell = [float(x) for x in cell[2:]]

    delta = .5
    lastImprovement = 0

    startDiff = None
    lastDiff = 9999

    i = -1
    while True:
        i += 1
        data = [float(x) for x in cell[2:]]
        jobs = [data]
        conDict = params[1]
        for p in params[0]:
            try:
                cons = conDict[p]
            except:
                cons = []
            j = data[:]
            j[p] -= delta
            for con in cons:
                j[con] = j[p]
            jobs.append(j)

            j = data[:]
            j[p] += delta
            for con in cons:
                j[con] = j[p]
            jobs.append(j)

        bestR = 999999
        bestRj = None

        bestmean = 999999
        bestmeanj = None

        bestW = lastDiff
        bestWj = 0
        bestRW = 9999

        numJobs = len(jobs)
        barLengths = 60
        for j, job in enumerate(jobs):
            progress = (j + 1) / numJobs
            progress = int(barLengths * progress)
            sys.stdout.write('\r Step {:3} ['.format(i + 1) + progress * '#' + (barLengths - progress) * '-' + ']')
            sys.stdout.flush()
            newCell = cell[:2] + ['{:7.4f}'.format(p) for p in job]  # + cell[5:]
            newCell = ' '.join(newCell) + '\n'
            reader['cell'] = newCell
            reader.write(fileName='work.ins')
            wR2, mean, weighted = evaluate('work')
            if wR2 < bestR:
                bestR = wR2
                bestRj = j
            if mean < bestmean:
                bestmean = mean
                bestmeanj = j
            if weighted < bestW:
                bestW = weighted
                bestWj = j
                bestRW = wR2

            if not startDiff:
                startDiff = weighted
        print()
        print()
        # jString = j2Name(bestWj)
        # print('  ', jString)
        print('   Old Cell:  ', cell2String(cell[2:], offset=15))
        print()
        # print('   Best wR2:            ', bestRj, cell2String(jobs[bestRj]))
        # print('   Best mean:           ', bestmeanj, cell2String(jobs[bestmeanj]))
        print('   New Cell:  ', cell2String(jobs[bestWj], offset=15))
        print()
        print('   DFIX Fit:     {:7.5f} / {:7.5f}\n'.format(bestW, startDiff))
        print('   Current wR2:  {}'.format('{:7.5f}\n'.format(bestRW) if bestRW < 10 else 'No imvprovements'))
        # input()
        cell = cell[:2] + ['{:7.4f}'.format(p) for p in jobs[bestWj]]  # + cell[5:]
        lastDiff = bestW
        if bestWj == 0:
            # print('No improvements found. Decreasing step size.')
            delta = delta / 2
            if delta < 0.005:
                print('Converged.')
                break
            if i - lastImprovement > 5:
                print('No improvements since 5 steps. Terminating.')
                break
        else:
            lastImprovement = i
        break
    print('\n\nOriginal Cell:', cell2String(originalCell, offset=15))
    print('   Final Cell:', cell2String(cell[2:], offset=15))

    print('\nOriginal DFIX fit: {:8.6f}'.format(startDiff))
    print('   Final DFIX fit: {:8.6f}'.format(bestW))


JDICT = {0: 'a',
         1: 'b',
         2: 'c',
         3: 'alpha',
         4: 'beta',
         5: 'gamma'}


# def j2Name(j):
#     if not j:
#         return 'No modifications.'
#     j -= 1
#     j = j // 2
#     name = ('Incremented ' if j % 2 else 'Decremented ') + JDICT[j]
#     return name


def cell2String(cell, offset=0):
    """
    Returns a nicely formatted string representation of unit cell parameters.
    :param cell: list
    :param offset: int<offset of row two and three>
    :return: str
    """
    cell = [float(x) for x in cell]
    return '{a:9.4f} {A:9.4f}\n{offset}{b:9.4f} {B:9.4f}\n{offset}{c:9.4f} {C:9.4f}'.format(a=cell[0],
                                                                                            b=cell[1],
                                                                                            c=cell[2],
                                                                                            A=cell[3],
                                                                                            B=cell[4],
                                                                                            C=cell[5],
                                                                                            offset=' ' * offset)



class Plotter(object):
    """
    Class for accumulating and plotting of data.
    """
    def __init__(self):

        self.values = {}

    def __call__(self, **kwargs):
        """
        Add data to the plotter.
        :param kwargs: Each key word will become a data series. Subsequent calls with equal keys will
         append the datum to the series.
        :return: None
        """
        for key, value in kwargs.items():
            try:
                self.values[key].append(value)
            except KeyError:
                self.values[key] = [value]

    def normalize(self):
        """
        Subtracts the first value of each series from each datum in the corresponding series.
        :return:
        """
        for key, values in self.values.items():
            m = values[0]
            self.values[key] = [v-m for v in values]


    def plot(self, x, y,  label=''):
        """
        Adds a data series to the plot.
        :param x: float
        :param y: float
        :param label: str
        :return: plot
        """
        plt.plot(x, y, marker='', label=label)
        plt.legend(loc='upper right')
        return plt

    def show(self):
        """
        Create and show a plot of the accumulated data.
        :return: None
        """
        if not PLOTAVAILABLE:
            print('Plot function not available. Please install matplotlib.'
                  'On some operating systems the python-tkinter package'
                  'is not installed by default. If installing matplotlib'
                  'does not solve this issue, pleas check if python-tkinter'
                  'is installed.')
            return
        self.normalize()
        plt = None
        for key, data in self.values.items():
            plt = self.plot(range(len(data)), data, label=key)
        plt.show()


class Array(object):
    def __init__(self, values):
        self.values = values

    def __iter__(self):
        for v in self.values:
            yield v

    def __len__(self):
        return len(self.values)

    def __add__(self, other):
        if type(other) == type(self):
            if not len(self) == len(other):
                raise ValueError('Arrays are not of equal length.')
            return Array([i + j for i, j in zip(self, other)])
        elif type(other) == float or type(other) == int:
            return Array([i + other for i in self.values])
        else:
            raise TypeError('Cannot add type Array to type {}.'.format(str(type(other))))

    def __imul__(self, other):
        if type(other) == int:
            self.values = [v * other for v in self.values]
            return self
        else:
            raise TypeError('Unsupported operation.')

    def __str__(self):
        return 'Array({})'.format(str(self.values))

    def __getitem__(self, val):
        try:
            start, stop, step = val.start, val.stop, val.step
        except AttributeError:
            pass
        else:
            return self.values[val]
        finally:
            return self.values[val]

    def dot(self, other):
        return sum([i * j for i, j in zip(self, other)])


class Matrix(object):
    def __init__(self, values):
        self.shape = (len(values[0]), len(values))
        self.values = values

    def __getitem__(self, val):
        if type(val) == tuple:
            if len(val) == 1:
                return self.values[val[0]]
            return self.values[val[1]][val[0]]

    def __str__(self):
        print(self.values)
        return '\n'.join([str(row) for row in self.values])

    def __imul__(self, other):
        if type(other) == int:
            self.values = [[v * other for v in row] for row in self.values]
            return self
        else:
            raise TypeError('Unsupported operation.')

    def transpose(self):
        rows = []
        for i in range(self.shape[0]):
            rows.append([r[i] for r in self.values])
        return Matrix(rows)

    def dot(self, other):
        newA = []
        for i, col in enumerate(self.transpose().values):
            s = sum([v * o for v, o in zip(col, other)])
            newA.append(s)
        return Array(newA)


class SymmetryElement(object):
    """
    Class representing a symmetry operation.
    """
    symm_ID = 1

    def __init__(self, symms, centric=False):
        """
        Constructor.
        """
        self.centric = centric
        self.symms = symms
        self.ID = SymmetryElement.symm_ID
        SymmetryElement.symm_ID += 1
        lines = []
        trans = []
        for symm in self.symms:
            line, t = self._parse_line(symm)
            lines.append(line)
            trans.append(t)
        self.matrix = Matrix(lines).transpose()
        self.trans = Array(trans)
        if centric:
            self.matrix *= -1
            self.trans *= -1

    def __str__(self):
        string = '''|{aa:2} {ab:2} {ac:2}|   |{v:2}|
|{ba:2} {bb:2} {bc:2}| + |{vv:2}|
|{ca:2} {cb:2} {cc:2}|   |{vvv:2}|'''.format(aa=self.matrix[0, 0],
                                             ab=self.matrix[0, 1],
                                             ac=self.matrix[0, 2],
                                             ba=self.matrix[1, 0],
                                             bb=self.matrix[1, 1],
                                             bc=self.matrix[1, 2],
                                             ca=self.matrix[2, 0],
                                             cb=self.matrix[2, 1],
                                             cc=self.matrix[2, 2],
                                             v=self.trans[0],
                                             vv=self.trans[1],
                                             vvv=self.trans[2])
        return string

    def __eq__(self, other):
        """
        Check two SymmetryElement instances for equivalence.
        Note that differences in lattice translation are ignored.
        :param other: SymmetryElement instance
        :return: True/False
        """
        m = (self.matrix == other.matrix).all()
        t1 = Array([v % 1 for v in self.trans])
        t2 = Array([v % 1 for v in other.trans])
        t = (t1 == t2).all()
        return m and t

    def __sub__(self, other):
        """
        Computes and returns the translational difference between two SymmetryElements. Returns 999.0 if the elements
        cannot be superimposed via an integer shift of the translational parts.
        :param other: SymmetryElement instance
        :return: float
        """
        if not self == other:
            return 999.
        return self.trans - other.trans

    def applyLattSymm(self, lattSymm):
        """
        Copies SymmetryElement instance and returns the copy after applying the translational part of 'lattSymm'.
        :param lattSymm: SymmetryElement.
        :return: SymmetryElement.
        """
        # newSymm = deepcopy(self)
        newSymm = SymmetryElement(self.toShelxl().split(','))
        newSymm.trans = Array([(self.trans[0] + lattSymm.trans[0]) / 1,
                         (self.trans[1] + lattSymm.trans[1]) / 1,
                         (self.trans[2] + lattSymm.trans[2]) / 1])
        newSymm.centric = self.centric
        return newSymm

    def toShelxl(self):
        """
        Generate and return string representation of Symmetry Operation in Shelxl syntax.
        :return: string.
        """
        axes = ['X', 'Y', 'Z']
        lines = []
        for i in range(3):
            op = self.matrix[i,]
            text = str(self.trans[i]) if self.trans[i] else ''
            for j in range(3):
                s = '' if not self.matrix[i, j] else axes[j]
                if self.matrix[i, j] < 0:
                    s = '-' + s
                elif s:
                    s = '+' + s
                text += s
            lines.append(text)
        return ', '.join(lines)

    def _parse_line(self, symm):
        symm = symm.lower().replace(' ', '')
        chars = ['x', 'y', 'z']
        line = []
        for char in chars:
            element, symm = self._partition(symm, char)
            line.append(element)
        if symm:
            trans = self._float(symm)
        else:
            trans = 0
        return line, trans

    def _float(self, string):
        try:
            return float(string)
        except ValueError:
            if '/' in string:
                string = string.replace('/', './') + '.'
                return eval('{}'.format(string))

    def _partition(self, symm, char):
        parts = symm.partition(char)
        if parts[1]:
            if parts[0]:
                sign = parts[0][-1]
            else:
                sign = '+'
            if sign is '-':
                return -1, ''.join((parts[0][:-1], parts[2]))
            else:
                return 1, ''.join((parts[0], parts[2])).replace('+', '')
        else:
            return 0, symm


class ShelxlLine(object):
    """
    Class representing a line in a Shelxl.res file.
    """

    def __init__(self, line, key=None):
        self.line = line
        self.key = key
        # print(self.line)

    def __str__(self):
        return 'LINE: ' + self.line

    def write(self):
        """
        Returns a string representation of a shelxl line as expected by SHELXL.
        :return: str
        """
        return self.line + '\n'


class ShelxlAtom(ShelxlLine):
    """
    Class Representing an Atom in a Shelxl.res file.
    """

    def __init__(self, line, virtual=False, key=None, resi=('', 0)):
        self.rawData = line
        self.key = key
        data = [word for word in line.split() if word]
        data = [float(word) if i else word for i, word in enumerate(data)]
        self.data = data
        self.resiClass = resi[1]
        if resi[1]:
            self.name = str(data[0]) + '_{}'.format(resi[0])
            # print(self.name)
        else:
            self.name = data[0]
        self.sfac = int(data[1])
        self.frac = Array(data[2:5])
        self.occ = (data[5] // 1, data[5] % 1)
        self.adp = Array(data[6:])
        if virtual:
            return
        if self.name[0].upper() == 'Q':
            self.qPeak = True
            ShelxlReader.CURRENTMOLECULE.addQPeak(self)
        else:
            self.qPeak = False
            ShelxlReader.CURRENTMOLECULE.addAtom(self)

    def __str__(self):
        return 'ATOM: {} {} {} {} {}'.format(self.name, self.sfac, self.frac, self.occ, self.adp)

    def write(self):
        """
        Returns a string representation of a shelxl atom as expected by SHELXL.
        :return: str
        """
        string = '{name:8} {sfac} {frac} {occ:6.3f} {adp}\n'.format(name=self.name.split('_')[0],
                                                                    sfac=self.sfac,
                                                                    frac=' '.join(
                                                                        ['{:6.4f}'.format(c) for c in self.frac]),
                                                                    occ=sum(self.occ),
                                                                    adp=' '.join(['{:6.4f}'.format(c) for c in self.adp]))
        if len(string) > 75:
            string = string.split()
            string = string[:7] + ['=\n   '] + string[7:] + ['\n']
            string = ' '.join(string)
        return string


class ShelxlMolecule(object):
    """
    Class representing a molecule-like object defined by the instructions given in a Shelxl.res file.
    Note that a class instance does not necessarily represent a chemical molecule. Applying symmetry operations might
    be necessary to create the whole chemical molecule, or mutiple chemical molecules can be represented by one
    class instance in cases with multiple molecules per asymmetric unit.
    """

    def __init__(self):
        self.sfacs = []
        self.customSfacData = {}
        self.atoms = []
        self.atomDict = OrderedDict()
        self.qPeaks = []
        self.cell = []
        self.cerr = []
        self.lattOps = []
        self.symms = []
        self.centric = False
        self.eqivs = {}
        self.dfixs = []
        self.dangs = []
        self.dfixErr = 0.02
        self.dangErr = 0.05
        self.eqivSymmMap = {}
        self.resiClass2Nums = {}
        self.resis = []

    def __iter__(self):
        for atom in self.atoms:
            yield atom

    def distance(self, atom1, atom2):
        """
        Compute the distance between two atoms with given names.
        Tanslational symmetry is ignored. The shortest distance between the given atoms is computed.
        :param atom1: str
        :param atom2: str
        :return: float
        """
        x, y, z = atom1.frac
        try:
            xx, yy, zz = atom2.frac + 99.5
        except TypeError:
            xx, yy, zz = Array(atom2.frac) + 99.5
        dx = (xx - x) % 1 - 0.5
        dy = (yy - y) % 1 - 0.5
        dz = (zz - z) % 1 - 0.5
        a, b, c, alpha, beta, gamma = self.cell
        alpha = alpha / 180. * pi
        beta = beta / 180. * pi
        gamma = gamma / 180. * pi
        dd = a ** 2 * dx ** 2 + b ** 2 * dy ** 2 + c ** 2 * dz ** 2 + 2 * b * c * cos(
            alpha) * dy * dz + 2 * a * c * cos(beta) * dx * dz + 2 * a * b * cos(gamma) * dx * dy
        return dd ** .5

    def asP1(self, full=False):
        """
        Generates and returns a new ShelxlMolecule instance where symmetry operations were applied to generate
        symmetry equivalent atoms. The returned instance represents the equivalent structural model in P1/P-1.
        :return: ShelxlMolecule
        """
        # full=True
        resiOffset = int(max(self.resis, key=lambda resi:int(resi[1]))[1])
        if resiOffset > 2:
            raise ValueError('Expanding structures to P1/P-1 is not supported for structures '
                             'containing multiple residues.')
        if not full:
            symms = [symm for symm in self.symms if not symm.centric]
        else:
            symms = self.symms[:]
        p1Atoms = {str(i + 2+resiOffset): [] for i in range(len(symms))}
        p1Mol = deepcopy(self)
        p1Mol.symms = []
        p1Mol.centric = False
        p1Mol.lattOps = []
        specials = []
        equivSymmMap = {}

        for atom in self.atoms:
            # p1Mol.addAtom(atom)
            for i, symm in enumerate(symms):
                resiKey = str(i + 2+resiOffset)
                newFrac = symm.matrix.dot(atom.frac)
                newFrac = newFrac + symm.trans
                vAtom = ShelxlAtom(atom.rawData, virtual=True)
                vAtom.frac = newFrac
                # vAtom.name += 'X{}'.format(i)
                vAtom.occ = (10, 1)
                atom.occ = (10, 1)
                distance = self.distance(atom, vAtom)
                specialName = atom.name + '_' + resiKey
                if distance < 0.1:
                    # print('Atom {} on special position. (SYMM {})'.format(atom.name, i))

                    specials.append(specialName)
                    continue
                else:
                    p1Atoms[resiKey].append(vAtom)
                    p1Mol.addAtom(vAtom)
                    p1Mol.atomDict[specialName] = vAtom

                # for key, eqiv in self.eqivs.items():
                #     if symm == eqiv:
                #         d = eqiv - symm
                #         newEqiv = SymmetryElement(['{}+X'.format(d[0]), '{}+Y'.format(d[1]), '{}+Z'.format(d[2])])
                #         if not key in equivSymmMap:
                #             equivSymmMap[key] = (resiKey, newEqiv)

        # Expand DFIX restraints.
        for atom in self.atoms:
            for i, symm in enumerate(symms):
                resiKey = str(i + 2+resiOffset)

                for k, dfix in enumerate(p1Mol.dfixs):
                    pairs = set(dfix[2])
                    newpairs = set()
                    for pair in pairs:
                        ii = 0
                        if '_$' in pair[0]:
                            ii = 0
                        elif '_$' in pair[1]:
                            ii = 1
                        if ii:
                            # eName = pair[ii]
                            # base, eId = eName.split('_')
                            # print(base, eId)
                            # try:
                            #     newID, newSymm = equivSymmMap[eId]
                            #     eName = base + '_{}_{}'.format(newID, eId)
                            #     print(eName)
                            # except KeyError:
                            #     pass
                            continue
                        if atom.name in pair:
                            newPair = (pair[0] + '_' + resiKey, pair[1] + '_' + resiKey)
                            if newPair[0] in specials or newPair[1] in specials:
                                continue
                            else:
                                newpairs.add(newPair)
                    p1Mol.dfixs[k] = (dfix[0], dfix[1], list(pairs.union(newpairs)))
            # input()
        # exit()
        self.eqivSymmMap = equivSymmMap
        # for eKey, newE in equivSymmMap.items():
        #     print(eKey)
        #     print(newE[0])
        #     print(newE[1])
        #     print()
        p1Mol._finalizeDfix()
        p1AtomList = []
        # exit()
        for key, values in p1Atoms.items():
            p1AtomList.append(ShelxlLine('RESI sym {}'.format(key)))
            # print( key)
            for value in values:
                p1AtomList.append(value)
                # print(value.name)
            # input()
        return p1Mol, p1AtomList

    def addAtom(self, atom):
        """
        Add an atom.
        :param atom: ShelxlAtom
        :return: None
        """
        self.atoms.append(atom)

    def addResidue(self, num, cls):
        try:
            self.resiClass2Nums[cls].append(num)
        except KeyError:
            self.resiClass2Nums[cls] = [num,]
        self.resis.append((cls,num))

    def addQPeak(self, qPeak):
        """
        Add a Q-Peak
        :param qPeak: ShelxlLine
        :return: None
        """
        self.qPeaks.append(qPeak)

    def setCell(self, cell):
        """
        Set the cell of the structural model.
        :param cell: list of six floats
        :return: None
        """
        self.cell = Array([float(c) for c in cell])

    def setWavelength(self, w):
        """
        Set the wavelengths.
        :param w: float
        :return: None
        """
        self.waveLength = float(w)

    def setZ(self, z):
        """
        Set Z.
        :param z: int
        :return: None
        """
        self.z = int(z)

    def setCerr(self, cerr):
        """
        Set the error of the unit cell parameters.
        :param cerr: list of six floats
        :return: None
        """
        self.cerr = Array([float(c) for c in cerr])

    def addSfac(self, sfac):
        """
        Add a predefined SFAC.
        :param sfac: str
        :return: None
        """
        self.sfacs.append(sfac)

    def addCustomSfac(self, data):
        """
        Add a custom SFAC
        :param data: list of strings
         eg.[ O,  0.237,  0.799,  0.649,  3.654,  0.800, 10.924,  0.261, 29.964,  0.035, 0.0, 0.0, 1.234, 0.660, 16.000]
        :return: None
        """
        symbol = data[0]
        data = [float(datum) for datum in data[1:]]
        data = [symbol] + data
        self.customSfacData[symbol] = data
        self.sfacs.append(symbol)

    def addSymm(self, symmData):
        """
        Add the content of a Shelxl SYMM command to generate the appropriate SymmetryElement instance.
        :param symmData: list of strings. eg.['1/2+X', '1/2+Y', '1/2+Z']
        :return: None
        """
        newSymm = SymmetryElement(symmData)
        self.symms.append(newSymm)
        for symm in self.lattOps:
            lattSymm = newSymm.applyLattSymm(symm)
            self.symms.append(lattSymm)
        if self.centric:
            self.symms.append(SymmetryElement(symmData, centric=True))
            for symm in self.lattOps:
                lattSymm = newSymm.applyLattSymm(symm)
                self.symms.append(lattSymm)

    def setCentric(self, value):
        """
        Defines the instance as representing a centrosymmetric structure. Generates the appropriate SymmetryElement
        instances automatically if called before adding further SYMM commands via self.addSymm().
        :param value: bool
        :return: None
        """
        self.centric = value
        self.symms.append(SymmetryElement(['-X', '-Y', '-Z']))
        self.symms[-1].centric = True

    def setLattOps(self, lattOps):
        """
        Adds lattice operations. If called before adding SYMM commands, the appropriate lattice operations are used
        automatically to generate further SymmetryElements.
        :param lattOps: list of SymmetryElement instances.
        :return: None
        """
        self.lattOps = lattOps

    def addDfix(self, value, err, atomPairs):
        """
        Adds a DFIX restraint to the model.
        :param value: float
        :param err: float or None
        :param atomPairs: list of tuples of strings
        :return: None
        """
        self.dfixs.append((value, err, [tuple([atomName.upper() for atomName in pair]) for pair in atomPairs]))

    def addDang(self, value, err, atomPairs):
        self.dangs.append((value, err, [tuple([atomName.upper() for atomName in pair]) for pair in atomPairs]))

    def addEqiv(self, name, data):
        """
        Adds an EQIV instruction to the model.
        :param name: str eg. '$1'
        :param data: list of strings equivalent to self.addSymm(symmData).
        :return: None
        """
        symm = SymmetryElement(data)
        self.eqivs[name] = symm

    def getAtom(self, atomName):
        """
        Returns ShelxlAtom instance with the given name
        :param atomName: string
        :return: ShelxlAtom instance
        """
        if '_$' in atomName:
            return self.getVirtualAtom(atomName)
        try:
            return self.atomDict[atomName]
        except KeyError:
            if '_' in atomName:
                base, cls = atomName.split('_')
                resiNums = self.resiClass2Nums[cls]
                names = ['{}_{}'.format(base, num) for num in resiNums]
                return [self.atomDict[name] for name in names if name in self.atomDict]
            # else:
            #     return [value for key, value in self.atomDict.items() if key.split('_')[0] == atomName]
            else:
                raise KeyError('No atom named {}.'.format(atomName))

    def getVirtualAtom(self, atomName):
        """
        Generates and returns a virtual atom -- an atom that is referenced with an EQIV instruction.
        :param atomName: str eg. 'O1_$1'
        :return: ShelxlAtom instance
        """
        base, equiv = atomName.split('_$')
        symm = self.eqivs['$' + equiv]
        atom = self.getAtom(base)
        newFrac = (symm.matrix.dot(atom.frac) + symm.trans)
        vAtom = ShelxlAtom(atom.rawData, virtual=True)
        vAtom.frac = newFrac
        return vAtom

    def finalize(self):
        """
        Called after reading a shelxl.res file. Sets up atom table and restraint table.
        :return: None
        """
        self.atomDict = OrderedDict()
        for atom in self.atoms:
            self.atomDict[atom.name] = atom
        self._finalizeDfix()
        # self.checkDfix()

    def checkDfix(self):
        """
        Compute and return the mean difference between restrained target values and actual distance as well as the
        weighted difference.
        :return: (float<mean>, float<weightedMean>)
        """
        if not self.dfixTable:
            raise ValueError('No DFIX restraints found.')
        vSum = 0
        wSum = 0
        sum = 0
        i = 0
        for atom1, dfixs in self.dfixTable.items():
            for atom2, data in dfixs.items():
                target, err = data
                a1s = self.getAtom(atom1)
                a2s = self.getAtom(atom2)
                if type(a1s) is list and type(a2s) is list:
                    for a1, a2 in zip(a1s, a2s):
                        d = self.distance(a1, a2)
                        diff = abs(d - target) ** 2
                        vSum += diff * err
                        wSum += err
                        sum += diff
                        i += 1
                elif type(a1s) is list or type(a2s) is list:
                    raise ValueError('Cellopt does not support restraints between different residues.')
                else:
                    d = self.distance(self.getAtom(atom1), self.getAtom(atom2))
                    diff = abs(d - target) ** 2
                    vSum += diff * err
                    wSum += err
                    sum += diff
                    i += 1
        # print((sum / i) ** .5, (vSum / wSum) ** .5)
        # exit()
        return (sum / i) ** .5, (vSum / wSum) ** .5

    def _finalizeDfix(self):
        dfixTable = {atom.name.upper(): {} for atom in self.atoms}
        for target, err, pairs in self.dfixs:
            if not err:
                err = self.dfixErr
            for atom1, atom2 in pairs:
                atom1 = atom1.upper()
                try:
                    tableRow1 = dfixTable[atom1]
                except KeyError:
                    tableRow1 = {}
                    dfixTable[atom1] = tableRow1
                try:
                    tableField1 = tableRow1[atom2]
                except KeyError:
                    tableRow1[atom2] = (target, err)

                atom2 = atom2.upper()
                try:
                    tableRow2 = dfixTable[atom2]
                except KeyError:
                    tableRow2 = {}
                    dfixTable[atom2] = tableRow2
                try:
                    tableField2 = tableRow2[atom1]
                except KeyError:
                    tableRow2[atom1] = (target, err)

        for target, err, pairs in self.dangs:
            if not err:
                err = self.dangErr
            for atom1, atom2 in pairs:
                atom1 = atom1.upper()
                try:
                    tableRow1 = dfixTable[atom1]
                except KeyError:
                    tableRow1 = {}
                    dfixTable[atom1] = tableRow1
                try:
                    tableField1 = tableRow1[atom2]
                except KeyError:
                    tableRow1[atom2] = (target, err)

                atom2 = atom2.upper()
                try:
                    tableRow2 = dfixTable[atom2]
                except KeyError:
                    tableRow2 = {}
                    dfixTable[atom2] = tableRow2
                try:
                    tableField2 = tableRow2[atom1]
                except KeyError:
                    tableRow2[atom1] = (target, err)
        self.dfixTable = dfixTable


class ShelxlReader(object):
    """
    Interface to read and interact with shelxl.res files.
    """
    CURRENTMOLECULE = None
    CURRENTINSTANCE = None

    def __init__(self):
        self.currentResi = ('', 0)
        self.lines = []
        self.atoms = []
        self._shelxlDict = {}

    def read(self, fileName):
        """
        Reads and parses a shelxl.res file with the given name. Returns the resulting ShelxlMolecule instance.
        :param fileName: str
        :return: ShelxlMolecule instance
        """
        ShelxlReader.CURRENTMOLECULE = ShelxlMolecule()
        ShelxlReader.CURRENTINSTANCE = self
        parser = LineParser()
        with Reader(fileName) as reader:
            for line in reader.readlines():
                if line[0] is '+':
                    reader.insert(line[1:-1])
                    line = '+    ' + line[1:]
                parser, line = parser(line)
                if line:
                    self.lines.append(line)
        # for line in self.lines:
        #     print(line)

        # for atom1 in self.CURRENTMOLECULE.atoms:
        # for atom2 in self.CURRENTMOLECULE.atoms:
        #     print(self.CURRENTMOLECULE.distance( atom1, atom2))
        # print(atom1.name)
        molecule = ShelxlReader.CURRENTMOLECULE
        molecule.finalize()
        ShelxlReader.CURRENTMOLECULE = None
        ShelxlReader.CURRENTINSTANCE = None
        self.molecule = molecule
        return molecule

    def write(self, fileName='out.res'):
        """
        Writes a shelxl.res file with the given filename representing the potentially modified structure.
        :param fileName: str
        :return: None
        """
        with open(fileName, 'w') as fp:
            for line in self.lines:
                key = line.key
                try:
                    data = self[key]
                except KeyError:
                    fp.write(line.write())
                else:
                    fp.write(data + '\n')

    def toP1(self, full=False):
        """
        Modifies the ShelxlReader instance to represent a structure expanded to P1/P-1. Centrosymmentric structures
        are expanded to P-1 by default. Override this behavior by setting full to True.
        :param full: bool
        :return: None
        """
        self.molecule, newAtoms = self.molecule.asP1(full=full)
        replacedDfix = False
        offset = 0
        for i, line in enumerate(self.lines):
            if line.key is 'latt':
                if '-' in line.line or full:
                    self.lines[i] = ShelxlLine('LATT -1')
                else:
                    self.lines[i] = ShelxlLine('LATT 1')
            if line.key is 'symm':
                self.lines[i] = ShelxlLine('')
            if line.key is 'hklf':
                self.lines = self.lines[:i + offset] + newAtoms + self.lines[i + offset:]
            if line.key is 'dfix':
                # self.lines[i] = ShelxlLine('')
                # print(line)
                if replacedDfix:
                    continue
                replacedDfix = True
                # print(len(self.molecule.dfixTable))
                dfixLines = []
                for value, err, allPairs in self.molecule.dfixs:
                    for l in range(len(allPairs) // 4):
                        l = 4 * l
                        pairs = allPairs[l:l + 4]
                        line = 'DFIX {value} {err} {pairs}'.format(value=value, err=err if err else '',
                                                                   pairs=' '.join(
                                                                       ['{} {}'.format(pair[0], pair[1]) for pair in
                                                                        pairs]))
                        dfixLines.append(ShelxlLine(line))
                        # print(line)
                before = len(self.lines)
                self.lines = self.lines[:i] + dfixLines + self.lines[i:]
                after = len(self.lines)
                offset = after - before
        # for line in self.lines:
        #     print(line.write())
        # self.write('p1.ins')
        # exit()

    def setCurrentResi(self, cls, num):
        """
        Sets the current RESIDUE.
        :param cls: str
        :param num: int
        :return: None
        """
        self.currentResi = (cls, num)

    def __getitem__(self, item):
        """
        Returns a structural attribute of the given name.
        :param item: str
        :return: object
        """
        return self._shelxlDict[item]

    def __setitem__(self, key, value):
        """
        Sets a structural attribute to a given value
        :param key: str
        :param value: object
        :return: None
        """
        self._shelxlDict[key] = value


class BaseParser(object):
    """
    Base class for parsers.
    """
    RETURNTYPE = None
    KEY = None

    def __init__(self, line):
        self.body = line

    def __call__(self, line):
        self.finished()
        return LineParser(), line

    def get(self, previousParser):
        if not self.body.endswith('='):
            self.finished()
            return previousParser, self.RETURNTYPE(self.body, key=self.KEY)
        else:
            self.body = self.body[:-1]
            return self, None

    def finished(self):
        pass


class LineParser(BaseParser):
    """
    Default parser for lines in shelxl.res files.
    The parser will identify if a more specialized parser is required, and creates one if necessary.
    """

    def __init__(self):
        self.COMMANDS = {'REM': self.doNothing,
                         'BEDE': self.doNothing,
                         'MOLE': self.doNothing,
                         'TITL': self.doNothing,
                         'CELL': CellParser,
                         'ZERR': CerrParser,
                         'SYMM': SymmParser,
                         'SFAC': SfacParser,
                         'UNIT': self.doNothing,
                         'TEMP': self.doNothing,
                         'L.S.': self.doNothing,
                         'BOND': self.doNothing,
                         'ACTA': self.doNothing,
                         'LIST': self.doNothing,
                         'PLAN': self.doNothing,
                         'WGHT': self.doNothing,
                         'FVAR': self.doNothing,
                         'SIMU': self.doNothing,
                         'RIGU': self.doNothing,
                         'SADI': self.doNothing,
                         'SAME': self.doNothing,
                         'DANG': DangParser,
                         'AFIX': self.doNothing,
                         'PART': self.doNothing,
                         'HKLF': HklfParser,
                         'ABIN': self.doNothing,
                         'ANIS': self.doNothing,
                         'ANSC': self.doNothing,
                         'ANSR': self.doNothing,
                         'BASF': self.doNothing,
                         'BIND': self.doNothing,
                         'BLOC': self.doNothing,
                         'BUMP': self.doNothing,
                         'CGLS': self.doNothing,
                         'CHIV': self.doNothing,
                         'CONF': self.doNothing,
                         'CONN': self.doNothing,
                         'DAMP': self.doNothing,
                         'DEFS': self.doNothing,
                         'DELU': self.doNothing,
                         'DFIX': DfixParser,
                         'DISP': self.doNothing,
                         'EADP': self.doNothing,
                         'EQIV': EqivParser,
                         'EXTI': self.doNothing,
                         'EXYZ': self.doNothing,
                         'FEND': self.doNothing,
                         'FLAT': self.doNothing,
                         'FMAP': self.doNothing,
                         'FRAG': self.doNothing,
                         'FREE': self.doNothing,
                         'GRID': self.doNothing,
                         'HFIX': self.doNothing,
                         'HTAB': self.doNothing,
                         'ISOR': self.doNothing,
                         'LATT': LattParser,
                         'LAUE': self.doNothing,
                         'MERG': self.doNothing,
                         'MORE': self.doNothing,
                         'MPLA': self.doNothing,
                         'NCSY': self.doNothing,
                         'NEUT': self.doNothing,
                         'OMIT': self.doNothing,
                         'PRIG': self.doNothing,
                         'RESI': ResiParser,
                         'RTAB': self.doNothing,
                         'SHEL': self.doNothing,
                         'SIZE': self.doNothing,
                         'SPEC': self.doNothing,
                         'STIR': self.doNothing,
                         'SUMP': self.doNothing,
                         'SWAT': self.doNothing,
                         'TWIN': self.doNothing,
                         'TWST': self.doNothing,
                         'WIGL': self.doNothing,
                         'WPDB': self.doNothing,
                         'XNPD': self.doNothing,
                         'Q': self.doNothing,
                         'END': self.doNothing,
                         'LONE': self.doNothing,
                         '+': self.doNothing,
                         }

    def __call__(self, line):
        line = line.rstrip('\n')
        if not line:
            return self.doNothing(line)
        command = line[:4]
        if command[0] is ' ':
            action = self.doNothing
        else:
            try:
                action = self.COMMANDS[command.rstrip()]
                if isinstance(action, type):
                    parser = action(line)
                    return parser.get(self)
            except KeyError:
                atomParser = AtomParser(line)
                return atomParser.get(self)
        return action(line)

    def doNothing(self, line):
        return self, ShelxlLine(line)


class AtomParser(BaseParser):
    """
    Parser of atom records in shelx.res files.
    """
    RETURNTYPE = ShelxlAtom

    def get(self, previousParser):
        if not self.body.endswith('='):
            self.finished()
            return previousParser, self.RETURNTYPE(self.body, key=self.KEY,
                                                   resi=ShelxlReader.CURRENTINSTANCE.currentResi)
        else:
            self.body = self.body[:-1]
            return self, None

    def __call__(self, line):
        return LineParser(), ShelxlAtom(self.body + line, resi=ShelxlReader.CURRENTINSTANCE.currentResi)


class CellParser(BaseParser):
    """
    Parser for CELL records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'cell'

    def finished(self):
        data = Array([float(word) for word in self.body.split()[1:] if word])
        ShelxlReader.CURRENTMOLECULE.setCell(data[1:])
        ShelxlReader.CURRENTMOLECULE.setWavelength(data[0])
        ShelxlReader.CURRENTINSTANCE['cell'] = self.body


class CerrParser(BaseParser):
    """
    Parser for CERR records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'cerr'

    def finished(self):
        data = Array([float(word) for word in self.body.split()[1:] if word])
        ShelxlReader.CURRENTMOLECULE.setCerr(data[1:])
        ShelxlReader.CURRENTMOLECULE.setZ(data[0])


class SfacParser(BaseParser):
    """
    Parser for SFAC records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine

    def __call__(self, line):

        self.body = self.body + '=\n' + line
        return LineParser(), ShelxlLine(self.body)

    def finished(self):
        custom = False
        words = [word for word in self.body.split()[1:] if word]
        for word in words:
            try:
                word = float(word)
            except ValueError:
                pass
            else:
                custom = True
                break
        if not custom:
            for sfac in words:
                ShelxlReader.CURRENTMOLECULE.addSfac(sfac)
        else:
            ShelxlReader.CURRENTMOLECULE.addCustomSfac(words)


class LattParser(BaseParser):
    """
    Parser for LATT records in shelxl.res files.
    """
    LATTDICT = {1: [],
                2: [SymmetryElement(('.5', '.5', '.5'))],
                3: [],
                4: [SymmetryElement(('.5', '.5', '0')),
                    SymmetryElement(('.5', '0', '.5')),
                    SymmetryElement(('0', '.5', '.5'))],
                5: [SymmetryElement(('0', '.5', '.5'))],
                6: [SymmetryElement(('.5', '0', '.5'))],
                7: [SymmetryElement(('.5', '.5', '0'))],
                }
    RETURNTYPE = ShelxlLine
    KEY = 'latt'

    def finished(self):
        data = [word for word in self.body.split() if word]
        latt = int(data[-1])
        if latt > 0:
            ShelxlReader.CURRENTMOLECULE.setCentric(True)
        lattOps = LattParser.LATTDICT[abs(latt)]
        ShelxlReader.CURRENTMOLECULE.setLattOps(lattOps)
        ShelxlReader.CURRENTINSTANCE['latt'] = self.body


class SymmParser(BaseParser):
    """
    Parser for SYMM records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'symm'

    def finished(self):
        symmData = self.body[4:].split(',')
        ShelxlReader.CURRENTMOLECULE.addSymm(symmData)


class DfixParser(BaseParser):
    """
    Parser for DFIX records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'dfix'

    def finished(self):
        data = [word for word in self.body.split() if word]
        command = data.pop(0)
        try:
            _, resi = command.split('_')
        except ValueError:
            resi = None
        value, data = float(data[0]), data[1:]
        try:
            err = float(data[0])
        except ValueError:
            err = None
        else:
            data = data[1:]
        pairs = []
        for i in range(len(data) // 2):
            i, j = 2 * i, 2 * i + 1
            name1, name2 = data[i], data[j]
            if resi:
                name1 = name1 + '_{}'.format(resi)
                name2 = name2 + '_{}'.format(resi)
            pairs.append((name1, name2))
            ShelxlReader.CURRENTMOLECULE.addDfix(value, err, pairs)


class DangParser(BaseParser):
    """
    Parser for DANG records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'dfix'

    def finished(self):
        data = [word for word in self.body.split() if word]
        command = data.pop(0)
        try:
            _, resi = command.split('_')
        except ValueError:
            resi = None
        value, data = float(data[0]), data[1:]
        try:
            err = float(data[0])
        except ValueError:
            err = None
        else:
            data = data[1:]
        pairs = []
        for i in range(len(data) // 2):
            i, j = 2 * i, 2 * i + 1
            name1, name2 = data[i], data[j]
            if resi:
                name1 = name1 + '_{}'.format(resi)
                name2 = name2 + '_{}'.format(resi)
            pairs.append((name1, name2))
            ShelxlReader.CURRENTMOLECULE.addDang(value, err, pairs)


class EqivParser(BaseParser):
    """
    Parser for EQIV records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine

    def finished(self):
        data = [word for word in self.body.split() if word][1:]
        name = data.pop(0)
        data = ' '.join(data)
        data = data.split(',')
        ShelxlReader.CURRENTMOLECULE.addEqiv(name, data)


class HklfParser(BaseParser):
    """
    Parser for HKLF records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'hklf'


class ResiParser(BaseParser):
    """
    Parser for RESI records in shelxl.res files.
    """
    RETURNTYPE = ShelxlLine
    KEY = 'resi'

    def finished(self):
        data = [word for word in self.body[4:].split() if word]
        cls, num = data[0], data[1]
        ShelxlReader.CURRENTINSTANCE.setCurrentResi(cls, num)
        ShelxlReader.CURRENTMOLECULE.addResidue(cls, num)


class Reader(object):
    """
    Super awesome class for reading files that might contain references to other files and you don't want to deal
    with that hassle.

    If file a.txt is:
        1
        2
        3
    and file b.txt is:
        a
        b
        c
    the code
        with Reader('a.txt') as reader:
            for line in reader.readlines():
                    if '2' in line:
                            reader.insert('b.txt')
                    if 'b' in line:
                            reader.remove()
                    print line
    will print
        1
        2
        a
        b
        3
    """

    def __init__(self, fileName):
        self.fileName = fileName
        self.inserted = None
        self.open = False
        self.fp = None

    def readlines(self, ):
        """
        Provides an interface equivalent to filePointer.readlines().
        :return: Yield string.
        """
        if not self.open:
            self.fp = open(self.fileName, 'r')
        while True:
            n = None
            if self.inserted:
                n = self.inserted.readline()
            if not n:
                n = self.fp.readline()
            if not n:
                raise StopIteration
            yield n

    def __exit__(self, *args):
        self.fp.close()
        try:
            self.inserted.close()
        except AttributeError:
            pass

    def __enter__(self):
        self.fp = open(self.fileName, 'r')
        return self

    def insert(self, fileName):
        """
        Insert a second file with a given name. After this method is called, each consecutive call to 'readlines' will
        yield a line of the inserted file until the inserted file yields EOF or 'remove' is called.
        :param fileName: string.
        :return: None
        """
        try:
            self.inserted = open(fileName, 'r')
        except FileNotFoundError:
            print('Cannot find insertion file.')

    def remove(self):
        """
        Removes a previously inserted file to stop yielding from the inserted file and continuing with the base file.
        :return:
        """
        self.inserted.close()
        self.inserted = None

    def fileInserted(self):
        """
        Check whether 'readlines' is currently yielding lines from an inserted file, or the base file.
        :return: bool.
        """
        return True if self.inserted else False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refine cell parameters against distance restraints.')
    parser.add_argument('fileName', type=str, help='Name of a shelxl result file.')
    parser.add_argument('-c', '--class', type=str, default=None,
                        help='Crystal class constraints for refinement.\nWARNING: The crystal class is ONLY used to'
                             'constrain the cell parameter refinement, and is NOT used to modify the structure'
                             'accordingly. Use this option with care.',
                        choices=['triclinic',
                                 'monoclinic',
                                 'orthorhombic',
                                 'tetragonal',
                                 'rhombohedral',
                                 'hexagonal',
                                 'cubic'])
    parser.add_argument('--expand', '-x', '-e', action='store_true',
                        help='Expand structure to P1 (P-1 for centric structures) before refinement. Note that a'
                             'potential center of invasion is not expanded to ensure the stability of intermediate'
                             'refinement steps. --class argument will be ignored.')
    parser.add_argument('--mode', '-m', type=str, default='default', nargs=1,
                        help="Specify optimization scheme. The {fast} scheme uses a simplex algorithm to optimize cell "
                             "parameters against DFIX restraints. The {default} scheme runs a SHELXL optimization step "
                             "after each time the simplex optimization converged and restarts the simplex (requires "
                             "SHELXL). The {accurate} scheme runs SHELXL as part of the simplex's evaluation step "
                             "(very slow, requires SHELXL)",
                        choices=['default', 'fast', 'accurate'])
    parser.add_argument('--plot', '-p', action='store_true',
                        help='Create diagnostic plot.')
    args = parser.parse_args()
    expand = args.expand
    crystalClass = args.__dict__['class']
    fileName = args.fileName

    if not os.path.isfile(fileName+'.res'):
        print('File {}.res is missing.'.format(fileName))
        exit(3)
    if not os.path.isfile(fileName+'.hkl'):
        print('File {}.hkl is missing.'.format(fileName))
        exit(4)
    mode = args.mode
    plot = args.plot
    if mode is 'default':
        run(fileName, p1=expand, overrideClass=crystalClass, plot=plot)
    elif mode is 'fast':
        run(fileName, p1=expand, overrideClass=crystalClass, fast=True, plot=plot)
    elif mode is 'accurate':
        run2(fileName, p1=expand)

    import urllib.request
    import json
    import subprocess
    try:
        localVersion = subprocess.check_output(["git2", "rev-parse", "HEAD"]).decode().strip()
    except:
        pass
    else:
        try:
            with urllib.request.urlopen('https://api.github.com/repos/JLuebben/CellOpt/commits/master ') as response:
                jsonData = response.read()
                data = json.loads(jsonData)
                remoteVersion = str(data['sha'])
                if not remoteVersion == localVersion:
                    print('A new version of cellopt.py is available at\n'
                          '     https://github.com/JLuebben/CellOpt')
        except:
            pass
