from cellopt.shelxlReader import ShelxlReader
from subprocess import call, STDOUT
from shutil import copyfile
import os
import sys
import argparse
from os.path import dirname, join

CLASSPARAMETERS = {'triclinic': ((0, 1, 2, 3, 4, 5), {}),
                   'monoclinic': ((0, 1, 2, 4), {}),
                   'orthorhombic': ((0, 1, 2), {}),
                   'tetragonal': ((0, 2), {1: 0}),
                   'rhombohedral': ((0, 3), {0: (1, 2), 3: (4, 5)}),
                   'hexagonal': ((0, 2), {0: (1,)}),
                   'cubic': ((0,), {0: (1, 2)})}


def callShelxl(fileName):
    FNULL = open(os.devnull, 'w')
    try:
        call(['shelxl.exe', fileName], stdout=FNULL, stderr=STDOUT)
    except:
        call(['shelxl', fileName], stdout=FNULL, stderr=STDOUT)


def evaluate(fileName):
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


def determineCrystalClass(cell):
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


def run(fileName, p1=False, overrideClass=None):
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

def j2Name(j):
    if not j:
        return 'No modifications.'
    j -= 1
    j = j//2
    name = ('Incremented ' if j%2 else 'Decremented ') + JDICT[j]
    return name


def cell2String(cell, offset=0):
    cell = [float(x) for x in cell]
    return '{a:9.4f} {A:9.4f}\n{offset}{b:9.4f} {B:9.4f}\n{offset}{c:9.4f} {C:9.4f}'.format(a=cell[0],
                                                                                            b=cell[1],
                                                                                            c=cell[2],
                                                                                            A=cell[3],
                                                                                            B=cell[4],
                                                                                            C=cell[5],
                                                                                            offset=' ' * offset)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refine cell parameters against distance restraints.')
    parser.add_argument('fileName', type=str, nargs=1, help='Name of a shelxl result file.')
    parser.add_argument('-c', '--class', type=str, default=None, nargs=1,
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
    parser.add_argument('--expand', action='store_true',
                        help='Expand structure to P1 before refinement. --class argument will be ignored.')
    args = parser.parse_args()
    expand = args.expand
    expand = True
    crystalClass = args.__dict__['class']
    fileName = args.fileName[0]
    # print(fileName)
    # crystalClass = 'monoclinic'
    run(fileName, p1=expand, overrideClass=crystalClass)
