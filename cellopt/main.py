from cellopt.shelxlReader import ShelxlReader
from subprocess import call, STDOUT
from shutil import copyfile
import os
from os.path import dirname,join


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
    molecule = reader.read(fileName+'.res')
    mean, weighted = molecule.checkDfix()
    return wR2, mean, weighted

def run(fileName):
    resFileName = fileName + '.res'
    fileDir = dirname(resFileName)
    copyfile(join(fileDir, fileName+'.hkl'), './work.hkl')
    reader = ShelxlReader()
    molecule = reader.read(resFileName)
    reader.toP1()
    # exit()
    # print(molecule.checkDfix())
    cell = reader['cell'].split()
    originalCell = [float(x) for x in cell[2:]]

    delta = .1
    lastImprovement = 0

    maxI = 100
    # for i in range(maxI):
    #     print('Step {}/{}'.format(i+1, maxI))
    i = -1
    while True:
        print('Step', i+1)
        i+=1
        a, b, c = [float(x) for x in cell[2:5]]
        params = (a, b, c)
        jobs = [params]
        for j, p in enumerate(params):
            job1 = list(params)
            job1[j] = job1[j] - delta
            jobs.append(job1)
            job2 = list(params)
            job2[j] = job2[j] + delta
            jobs.append(job2)

        bestR = 999999
        bestRj = None

        bestmean = 999999
        bestmeanj = None

        bestW = 999999
        bestWj = None
        for j, job in enumerate(jobs):
            print(' substep: {}/{}'.format( j+1, len(jobs)))
            newCell = cell[:2] + ['{:7.4f}'.format(p) for p in job] + cell[5:]
            newCell = ' '.join(newCell) + '\n'
            # print(newCell)
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

        print('   Old Cell:              ', [float(x) for x in cell[2:5]])
        print('   Best wR2:            ', bestRj, jobs[bestRj])
        print('   Best mean:           ', bestmeanj, jobs[bestmeanj])
        print('  !Best weighted mean:  ',  bestWj, jobs[bestWj], '!')
        # input()
        cell = cell[:2] + ['{:7.4f}'.format(p) for p in jobs[bestWj]] + cell[5:]
        if bestWj == 0:
            print('No improvements found. Decreasing step size.')
            delta = delta/2
            if delta < 0.005:
                print('Converged.')
                break
            if i-lastImprovement>5:
                print('No improvements since 5 steps. Terminating.')
                break
        else:
            lastImprovement = i
        # print(reader['cell'])
    print('\n\nOriginal Cell:', originalCell)
    print('   Final Cell:', [float(x) for x in cell[2:]])


if __name__ == '__main__':
    from sys import argv
    run('s1')
