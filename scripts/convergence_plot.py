import subprocess
import collections
import sys

Error = collections.namedtuple('Error', 'dx L_inf L_1 L_2')


def parseError(line):
    tokens = line.split()
    dx = float(tokens[tokens.index('dx:')+1])
    L_inf = float(tokens[tokens.index('L_inf:')+1])
    L_1 = float(tokens[tokens.index('L_1:')+1])
    L_2 = float(tokens[tokens.index('L_2:')+1])
    return Error(dx=dx, L_inf=L_inf, L_1=L_1, L_2=L_2)


def toCSV(errors, fname):
    with open(fname, "w") as f:
        print("dx,L_inf,L_1,L_2", file=f)
        for error in errors:
            print("%f,%f,%f,%f" %
                  (error.dx, error.L_inf, error.L_1, error.L_2), file=f)


divErrors = []
rotErrors = []
lapErrors = []

numArg = len(sys.argv)
if (not numArg == 2):
    print("intended usage %s <path/to/stenciLBinary>" % (sys.argv[0]))
    sys.exit(1)

processName = sys.argv[1]

for ny in range(16, 128+1, 16):
    errors = subprocess.run([processName, str(ny)],
                            stdout=subprocess.PIPE).stdout.decode('utf-8')

    for line in errors.splitlines():
        if (line.startswith('[div]')):
            divErrors.append(parseError(line))
        elif (line.startswith('[rot]')):
            rotErrors.append(parseError(line))
        elif (line.startswith('[lap]')):
            lapErrors.append(parseError(line))
        elif (line.startswith('---')):
            continue
        else:
            print("unexpected line in ouptut")
            sys.exit(1)

toCSV(divErrors, processName[2:] + "ConvDiv.csv")
toCSV(divErrors, processName[2:] + "ConvRot.csv")
toCSV(divErrors, processName[2:] + "ConvLap.csv")
