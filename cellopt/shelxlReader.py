import numpy as np
from copy import deepcopy


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
        self.matrix = np.matrix(lines).transpose()
        self.trans = np.array(trans)
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
        m = (self.matrix == other.matrix).all()
        t1 = np.array([v % 1 for v in self.trans])
        t2 = np.array([v % 1 for v in other.trans])
        t = (t1 == t2).all()
        return m and t

    def applyLattSymm(self, lattSymm):
        """
        Copies SymmetryElement instance and returns the copy after applying the translational part of 'lattSymm'.
        :param lattSymm: SymmetryElement.
        :return: SymmetryElement.
        """
        # newSymm = deepcopy(self)
        newSymm = SymmetryElement(self.toShelxl().split(','))
        newSymm.trans = [(self.trans[0] + lattSymm.trans[0]) / 1,
                         (self.trans[1] + lattSymm.trans[1]) / 1,
                         (self.trans[2] + lattSymm.trans[2]) / 1]
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
    def __init__(self, line, key=None):
        self.line = line
        self.key = key
        # print(self.line)

    def __str__(self):
        return 'LINE: ' + self.line

    def write(self):
        return self.line + '\n'


class ShelxlAtom(object):
    def __init__(self, line, virtual=False, key=None, resi=('', 0)):
        self.rawData = line
        self.key = key
        data = [word for word in line.split() if word]
        data = [float(word) if i else word for i, word in enumerate(data)]
        self.data = data
        self.resiClass = resi[0]
        # print(resi)
        if resi[1]:
            self.name = data[0]+'_{}'.format(resi[1])
            # print(self.name)
        else:
            self.name = data[0]
        self.sfac = int(data[1])
        self.frac = np.array(data[2:5])
        self.occ = (data[5] // 1, data[5] % 1)
        self.adp = np.array(data[6:])
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
        return '{name:8} {sfac} {frac} {occ:6.3f} {adp}\n'.format(name=self.name,
                                                                  sfac=self.sfac,
                                                                  frac=' '.join(
                                                                      ['{:6.4f}'.format(c) for c in self.frac]),
                                                                  occ=sum(self.occ),
                                                                  adp=' '.join(['{:6.4f}'.format(c) for c in self.adp]))


class ShelxlMolecule(object):
    def __init__(self):
        self.sfacs = []
        self.customSfacData = {}
        self.atoms = []
        self.atomDict = {}
        self.qPeaks = []
        self.cell = []
        self.cerr = []
        self.lattOps = []
        self.symms = []
        self.centric = False
        self.eqivs = {}
        self.dfixs = []
        self.dfixErr = 0.02

    def __iter__(self):
        for atom in self.atoms:
            yield atom

    def distance(self, atom1, atom2):
        x, y, z = atom1.frac
        try:
            xx, yy, zz = atom2.frac + 99.5
        except TypeError:
            xx, yy, zz = np.array(atom2.frac) + 99.5
        dx = (xx - x) % 1 - 0.5
        dy = (yy - y) % 1 - 0.5
        dz = (zz - z) % 1 - 0.5
        a, b, c, alpha, beta, gamma = self.cell
        alpha = alpha / 180. * np.pi
        beta = beta / 180. * np.pi
        gamma = gamma / 180. * np.pi
        dd = a ** 2 * dx ** 2 + b ** 2 * dy ** 2 + c ** 2 * dz ** 2 + 2 * b * c * np.cos(
            alpha) * dy * dz + 2 * a * c * np.cos(beta) * dx * dz + 2 * a * b * np.cos(gamma) * dx * dy
        return dd ** .5

    def asP1(self):
        symms = [symm for symm in self.symms if not symm.centric]
        p1Atoms = {i + 2: [] for i in range(len(symms))}
        p1Mol = deepcopy(self)
        p1Mol.symms = []
        p1Mol.centric = False
        p1Mol.lattOps = []
        specials = []

        for atom in self.atoms:
            # p1Mol.addAtom(atom)
            for i, symm in enumerate(symms):
                resiKey = str(i + 2)
                newFrac = (np.dot(atom.frac, symm.matrix) + symm.trans).flatten().tolist()[0]
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
                else:
                    p1Atoms[i + 2].append(vAtom)
                    p1Mol.addAtom(vAtom)
                    p1Mol.atomDict[specialName] = vAtom

        #         for key, eqiv in self.eqivs.items():
        #             if symm == eqiv:
        #                 print(i, key)
        #                 print(symm)
        #                 print(eqiv)
        #                 print()
        #         print()

                for k, dfix in enumerate(p1Mol.dfixs):
                    pairs = set(dfix[2])
                    # print(dfix)
                    newpairs = set()
                    for pair in pairs:
                        if atom.name in pair:
                            # print(atom.name, pair)
                            newPair = (pair[0]+'_'+resiKey, pair[1]+'_'+resiKey)
                            if newPair[0] in specials or newPair[1] in specials:
                                continue
                            newpairs.add(newPair)
                    p1Mol.dfixs[k] = (dfix[0], dfix[1], list(pairs.union(newpairs)))
        p1Mol._finalizeDfix()
        p1AtomList = []
        for key, values in p1Atoms.items():
            p1AtomList.append(ShelxlLine('RESI sym {}'.format(key)))
            for value in values:
                p1AtomList.append(value)
        return p1Mol, p1AtomList

    def addAtom(self, atom):
        self.atoms.append(atom)

    def addQPeak(self, qPeak):
        self.qPeaks.append(qPeak)

    def setCell(self, cell):
        self.cell = np.array([float(c) for c in cell])

    def setWavelength(self, w):
        self.waveLength = float(w)

    def setZ(self, z):
        self.z = int(z)

    def setCerr(self, cerr):
        self.cerr = np.array([float(c) for c in cerr])

    def addSfac(self, sfac):
        self.sfacs.append(sfac)

    def addCustomSfac(self, data):
        symbol = data[0]
        data = [float(datum) for datum in data[1:]]
        data = [symbol] + data
        self.customSfacData[symbol] = data
        self.sfacs.append(symbol)

    def addSymm(self, symmData):
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
        self.centric = value
        self.symms.append(SymmetryElement(['-X', '-Y', '-Z']))
        self.symms[-1].centric=True

    def setLattOps(self, lattOps):
        self.lattOps = lattOps

    def addDfix(self, value, err, targets):
        self.dfixs.append((value, err, targets))

    def addEqiv(self, name, data):
        symm = SymmetryElement(data)
        self.eqivs[name] = symm

    def getAtom(self, atomName):
        if '_$' in atomName:
            return self.getVirtualAtom(atomName)
        return self.atomDict[atomName]

    def getVirtualAtom(self, atomName):
        base, equiv = atomName.split('_$')
        symm = self.eqivs['$' + equiv]
        atom = self.getAtom(base)
        newFrac = (np.dot(atom.frac, symm.matrix) + symm.trans).flatten().tolist()[0]
        vAtom = ShelxlAtom(atom.rawData, virtual=True)
        vAtom.frac = newFrac
        return vAtom

    def finalize(self):
        self.atomDict = {}
        for atom in self.atoms:
            self.atomDict[atom.name] = atom
        self._finalizeDfix()
        # self.checkDfix()

    def checkDfix(self):
        vSum = 0
        wSum = 0
        sum = 0
        i = 0
        for atom1, dfixs in self.dfixTable.items():
            for atom2, data in dfixs.items():
                target, err = data
                try:
                    d = self.distance(self.getAtom(atom1), self.getAtom(atom2))
                except ValueError:
                    print(atom1, atom2)
                diff = abs(d - target) ** 2
                vSum += diff * err
                wSum += err
                sum += diff
                i += 1

        return (sum / i) ** .5, (vSum / wSum) ** .5
        # print('         Mean squared difference: {:6.4f}'.format((sum/i)**.5))
        # print('Weighted mean squared difference: {:6.4f}'.format((vSum/wSum)**.5))

    def _finalizeDfix(self):
        dfixTable = {atom.name.upper(): {} for atom in self.atoms}
        # for atom in self.atoms:
        #     print(atom.name)
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
                # else:
                #     tableField1.append((target, err))

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
                # else:
                #     tableField2.append((target, err))
        self.dfixTable = dfixTable


class ShelxlReader(object):
    CURRENTMOLECULE = None
    CURRENTINSTANCE = None

    def __init__(self):
        self.currentResi = ('', 0)
        self.lines = []
        self.atoms = []
        self._shelxlDict = {}

    def read(self, fileName):
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
        with open(fileName, 'w') as fp:
            for line in self.lines:
                key = line.key
                try:
                    data = self[key]
                except KeyError:
                    fp.write(line.write())
                else:
                    fp.write(data + '\n')

    def toP1(self):
        self.molecule, newAtoms = self.molecule.asP1()
        replacedDfix = False
        offset = 0
        for i, line in enumerate(self.lines):
            if line.key is 'latt':
                if '-' in line.line:
                    self.lines[i] = ShelxlLine('LATT -1')
                else:
                    self.lines[i] = ShelxlLine('LATT 1')
            if line.key is 'symm':
                self.lines[i] = ShelxlLine('')
            if line.key is 'hklf':
                self.lines = self.lines[:i+offset] + newAtoms + self.lines[i+offset:]
            if line.key is 'dfix':
                # self.lines[i] = ShelxlLine('')
                # print(line)
                if replacedDfix:

                    continue
                replacedDfix = True
                # print(len(self.molecule.dfixTable))
                dfixLines = []
                for value, err, allPairs in self.molecule.dfixs:
                    for l in range(len(allPairs)//4):
                        l = 4*l
                        pairs = allPairs[l:l+4]
                        line = 'DFIX {value} {err} {pairs}'.format(value=value, err=err if err else '',
                                                                   pairs=' '.join(['{} {}'.format(pair[0], pair[1]) for pair in pairs]))
                        dfixLines.append(ShelxlLine(line))
                        # print(line)
                before = len(self.lines)
                self.lines = self.lines[:i] + dfixLines + self.lines[i:]
                after = len(self.lines)
                offset = after-before
        # for line in self.lines:
        #     print(line.write())
        # self.write('p1.ins')
        # exit()

    def setCurrentResi(self, cls, num):
        self.currentResi = (cls, num)

    def __getitem__(self, item):
        return self._shelxlDict[item]

    def __setitem__(self, key, value):
        self._shelxlDict[key] = value


class BaseParser(object):
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
    def __init__(self):
        self.COMMANDS = {'REM': self.doNothing,
                         'BEDE': self.doNothing,
                         'MOLE': self.doNothing,
                         'TITL': self.doNothing,
                         'CELL': CellParser,
                         'ZERR': CerrParser,
                         'LATT': self.doNothing,
                         'SYMM': SymmParser,
                         'SFAC': SfacParser,
                         'UNIT': self.doNothing,
                         'TEMP': self.doNothing,
                         'L.S.': self.doNothing,
                         'BOND': self.doNothing,
                         'ACTA': self.doNothing,
                         'LIST': self.doNothing,
                         'FMAP': self.doNothing,
                         'PLAN': self.doNothing,
                         'WGHT': self.doNothing,
                         'FVAR': self.doNothing,
                         'SIMU': self.doNothing,
                         'RIGU': self.doNothing,
                         'SADI': self.doNothing,
                         'SAME': self.doNothing,
                         'DANG': self.doNothing,
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
                         'PLAN': self.doNothing,
                         'PRIG': self.doNothing,
                         'RESI': ResiParser,
                         'RTAB': self.doNothing,
                         'SADI': self.doNothing,
                         'SAME': self.doNothing,
                         'SHEL': self.doNothing,
                         'SIMU': self.doNothing,
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
                         'REM': self.doNothing,
                         'Q': self.doNothing,
                         'END': self.doNothing,
                         'BEDE': self.doNothing,
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
    RETURNTYPE = ShelxlAtom

    def get(self, previousParser):
        if not self.body.endswith('='):
            self.finished()
            return previousParser, self.RETURNTYPE(self.body, key=self.KEY, resi=ShelxlReader.CURRENTINSTANCE.currentResi)
        else:
            self.body = self.body[:-1]
            return self, None

    def __call__(self, line):
        return LineParser(), ShelxlAtom(self.body + line, resi=ShelxlReader.CURRENTINSTANCE.currentResi)


class CellParser(BaseParser):
    RETURNTYPE = ShelxlLine
    KEY = 'cell'

    def finished(self):
        data = np.array([float(word) for word in self.body.split()[1:] if word])
        ShelxlReader.CURRENTMOLECULE.setCell(data[1:])
        ShelxlReader.CURRENTMOLECULE.setWavelength(data[0])
        ShelxlReader.CURRENTINSTANCE['cell'] = self.body


class CerrParser(BaseParser):
    RETURNTYPE = ShelxlLine
    KEY = 'cerr'

    def finished(self):
        data = np.array([float(word) for word in self.body.split()[1:] if word])
        ShelxlReader.CURRENTMOLECULE.setCerr(data[1:])
        ShelxlReader.CURRENTMOLECULE.setZ(data[0])


class SfacParser(BaseParser):
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
    RETURNTYPE = ShelxlLine
    KEY = 'symm'

    def finished(self):
        symmData = self.body[4:].split(',')
        ShelxlReader.CURRENTMOLECULE.addSymm(symmData)


class DfixParser(BaseParser):
    RETURNTYPE = ShelxlLine
    KEY = 'dfix'

    def finished(self):
        data = [word for word in self.body[4:].split() if word]
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
            pairs.append((data[i], data[j]))
            ShelxlReader.CURRENTMOLECULE.addDfix(value, err, pairs)


class EqivParser(BaseParser):
    RETURNTYPE = ShelxlLine

    def finished(self):
        data = [word for word in self.body.split() if word][1:]
        name = data.pop(0)
        data = ' '.join(data)
        data = data.split(',')
        ShelxlReader.CURRENTMOLECULE.addEqiv(name, data)


class HklfParser(BaseParser):
    RETURNTYPE = ShelxlLine
    KEY = 'hklf'


class ResiParser(BaseParser):
    RETURNTYPE = ShelxlLine
    KEY = 'resi'

    def finished(self):
        data = [word for word in self.body[4:].split() if word]
        cls, num = data[0], data[1]
        ShelxlReader.CURRENTINSTANCE.setCurrentResi(cls, num)


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
    shelxlFile = ShelxlReader()
    molecule = shelxlFile.read('s1.ins')
    print(molecule.checkDfix())
