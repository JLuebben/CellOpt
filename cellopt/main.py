from cellopt.shelxlReader import ShelxlReader
from subprocess import call, STDOUT
from shutil import copyfile
import os
import sys
import argparse
from os.path import dirname, join
try:
    import matplotlib.pyplot as plt
except ImportError:
    PLOTAVAILABLE = False
else:
    PLOTAVAILABLE = True

CLASSPARAMETERS = {'triclinic': ((0, 1, 2, 3, 4, 5), {}),
                   'monoclinic': ((0, 1, 2, 4), {}),
                   'orthorhombic': ((0, 1, 2), {}),
                   'tetragonal': ((0, 2), {1: 0}),
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
    mean, weighted = molecule.checkDfix()
    return wR2, mean, weighted


def quickEvaluate(molecule, cell):
    """
    Evaluates the DFIX fit to a given cell
    :param molecule: ShelxlMolecule instance
    :param cell: list of six floats
    :return: float<meanDfixFit>, float<weightedDfixFit>
    """
    molecule.cell = cell
    return molecule.checkDfix()


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
        reader.toP1()
        cls = 'triclinic'
        params = CLASSPARAMETERS[cls]
    originalCell = [float(x) for x in cell[2:]]
    startDiff, _ = molecule.checkDfix()
    startDiff0 =startDiff
    lastDiff = 9999

    iterations = 25

    i = -1
    barLengths = 30

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
            progress = (i) / iterations
            progress = int(barLengths * progress)
            sys.stdout.write(
                '\r [' + progress * '#' + (barLengths - progress) * '-' + '] {fit:8.6f} {cell}'.format(fit=weighted,
                                                                                                       cell=' '.join([
                                                                                                           '{:9.4f}'.format(
                                                                                                               p)
                                                                                                           for
                                                                                                           p
                                                                                                           in
                                                                                                           job])))
            sys.stdout.flush()
        else:
            progress = barLengths
            sys.stdout.write(
                '\r [' + progress * '#' + (barLengths - progress) * '-' + '] {fit:8.6f} {cell}'.format(fit=weighted,
                                                                                                       cell=' '.join([
                                                                                                           '{:9.4f}'.format(
                                                                                                               p)
                                                                                                           for
                                                                                                           p
                                                                                                           in
                                                                                                           job])))
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
            # m = max(values)
            m = values[0]
            # self.values[key] = [v/m for v in values]
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
        self.normalize()
        plt = None
        for key, data in self.values.items():
            plt = self.plot(range(len(data)), data, label=key)
        plt.show()


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
    # expand = True
    crystalClass = args.__dict__['class']
    # crystalClass = 'orthorhombic'
    fileName = args.fileName
    # print(fileName)
    # crystalClass = 'monoclinic'
    mode = args.mode
    # mode = 'fast'
    # mode = 'accurate'
    plot = args.plot
    if mode is 'default':
        run(fileName, p1=expand, overrideClass=crystalClass, plot=plot)
    elif mode is 'fast':
        run(fileName, p1=expand, overrideClass=crystalClass, fast=True, plot=plot)
    elif mode is 'accurate':
        run2(fileName, p1=expand)
