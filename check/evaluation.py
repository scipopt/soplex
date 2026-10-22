#! /usr/bin/env python3

import sys
import os
import re
import json
import fractions
import io

# this function checks for the type of parsed string
def typeofvalue(text):
    try:
        int(text)
        return int
    except ValueError:
        pass

    try:
        float(text)
        return float
    except ValueError:
        pass

    try:
        fractions.Fraction(text)
        return fractions.Fraction
    except ValueError:
        pass

    return str

def isnum(text):
    return typeofvalue(text) in [int, float, fractions.Fraction]

def fracttofloat(value):
    try:
        return float(value)
    except ValueError:
        return float(fractions.Fraction(value))

# convert a numeric token to float, keep non-numerics (e.g. '-') as is
def tofloat(value):
    return fracttofloat(value) if isnum(value) else value

def fieldval(line):
    return line.split(':', 1)[1].split()[0]

# the larger of two values, ignoring non-numeric entries
def vmax(a, b):
    an = isinstance(a, (int, float))
    bn = isinstance(b, (int, float))
    if not an and not bn:
        return '-'
    if not an:
        return b
    if not bn:
        return a
    return a if a > b else b

if len(sys.argv) > 2:
    if not len(sys.argv) == 3:
        print('usage (for decomposition solve output): '+sys.argv[0]+'-d <soplex_test_run>.out')
        quit()
elif not len(sys.argv) == 2:
    print('usage: '+sys.argv[0]+' <soplex_test_run>.out')
    quit()

# if the -d flag is provided, then we must produce the decomposition solve output
decomp = False
if sys.argv[1] == '-d':
    decomp = True

decompcolumns = ['ppiv','dpiv','avgpdegen','avgddegen','redprob','compprob','rediter','compiter','algiter']
perfcolumns = ['instructions','cycles','clock','nmsfreq','normtime']

if decomp:
    outname = sys.argv[2]
else:
    outname = sys.argv[1]
dataname = outname.replace('.out','.json')

outfile = io.open(outname,'r',encoding='utf-8',errors='ignore')
outlines = outfile.readlines()
outfile.close()

testset = outname.split('/')[-1].split('.')[1]
settings = outname.split('/')[-1].split('.')[-2]

# locate test file
testname = ''
for d in ['testset', 'instancedata/testsets']:
    if os.path.isfile(d+'/'+testset+'.test'):
        testname = d+'/'+testset+'.test'
        break

# locate solution file
soluname = ''
for f in [testset, testset.split('_', 1)[0], testset.split('-', 1)[0], 'all']:
    for d in ['testset', 'instancedata/testsets']:
        if os.path.isfile(d+'/'+f+'.solu'):
            soluname = d+'/'+f+'.solu'
            break
    if soluname:
        break

# maximum length of instance names
namelength = 18

# tolerances for solution value check
tolerance = 1e-6
reltol = 1e-5

# solver statuses indicating reached limit
limitstats = ('timeout', 'iterlimit', 'objlimit')

instances = {}
objsenses = {}
instancename = ''
section = 'soplex'
stats = False
failflag = False
inorigcol = False
printperf = False
printedIdentifier = False

for outline in outlines:
    # print identifier
    if outline.startswith('@01'):
        # convert line to used instance name
        linesplit = outline.split('/')
        linesplit = linesplit[-1].rstrip(' \n').rstrip('.gz').rstrip('.GZ').rstrip('.z').rstrip('.Z')
        linesplit = linesplit.split('.')
        instancename = linesplit[0]
        for i in range(1, len(linesplit)-1):
            instancename = instancename + '.' + linesplit[i]
        length = len(instancename)
        if length > namelength:
            shortname = instancename[0:int(namelength/2)-1] + '~' + instancename[length-int(namelength/2):]
        else:
            shortname = instancename

        # initialize new data set
        instances[instancename] = {}
        instances[instancename]['status'] = 'fail (abort)'
        instances[instancename]['name'] = shortname
        instances[instancename]['settings'] = settings
        instances[instancename]['testset'] = testset
        # wait for statistics block
        section = 'soplex'
        stats = False
        failflag = False
        inorigcol = False
        continue

    # invalidate instancename
    if outline.startswith('=ready='):
        section = 'soplex'
        instancename = ''
        stats = False
        continue

    if outline.startswith('=perplex='):
        section = 'perplex'
        stats = False
        instances[instancename]['perplex'] = 'unknown'
        continue

    if outline.startswith('=qsoptex='):
        section = 'qsoptex'
        stats = False
        instances[instancename]['qso:stat'] = 'unknown'
        instances[instancename]['qso:lpval'] = '--'
        instances[instancename]['qso:time'] = instances[instancename]['timelimit']
        instances[instancename]['qso:prec'] = 64
        continue

    if instancename == '':
        continue

    # frequency normalization markers injected by evalcheck_cluster.sh
    if outline.startswith('@06'):
        for token in outline.split()[1:]:
            key, _, val = token.partition('=')
            if key == 'instructions:u' and isnum(val):
                instances[instancename]['instructions'] = int(fracttofloat(val))
            elif key == 'cpu-cycles:u' and isnum(val):
                instances[instancename]['cycles'] = int(fracttofloat(val))
            elif key == 'task-clock' and isnum(val):
                instances[instancename]['clock'] = fracttofloat(val)
        printperf = True
        continue

    if outline.startswith('@10'):
        for token in outline.split()[1:]:
            key, _, val = token.partition('=')
            if key == 'nmsfreq' and isnum(val):
                instances[instancename]['nmsfreq'] = int(fracttofloat(val))
            elif key == 'normtime' and isnum(val):
                instances[instancename]['normtime'] = fracttofloat(val)
        printperf = True
        continue

    if section == 'perplex':
        if outline.find('No such file or directory') >= 0:
            instances[instancename]['perplex'] = 'readerror'
        elif outline.startswith('Basis read'):
            instances[instancename]['perplex'] = 'timeout'
        elif outline.startswith('Solution is optimal'):
            if instances[instancename]['perplex'] == 'dinfeas':
                instances[instancename]['perplex'] = 'optimal'
            elif instances[instancename]['perplex'] == 'pdinfeas' or instances[instancename]['perplex'] == 'timeout':
                instances[instancename]['perplex'] = 'pinfeas'
        elif outline.startswith('Solution is feasible'):
            if instances[instancename]['perplex'] == 'pinfeas':
                instances[instancename]['perplex'] = 'optimal'
            elif instances[instancename]['perplex'] == 'pdinfeas' or instances[instancename]['perplex'] == 'timeout':
                instances[instancename]['perplex'] = 'dinfeas'
        elif outline.startswith('Solution is not ') and instances[instancename]['perplex'] == 'timeout':
            instances[instancename]['perplex'] = 'pdinfeas'
        continue

    if section == 'qsoptex':
        if outline.startswith('Time for SOLVER:'):
            instances[instancename]['qso:time'] = fracttofloat(outline.split()[-2])
        elif outline.find('Problem Solved Exactly') >= 0:
            instances[instancename]['qso:stat'] = 'optimal'
        elif outline.find('Problem Is Infeasible') >= 0:
            instances[instancename]['qso:stat'] = 'infeasible'
        elif outline.startswith('@24') and instances[instancename]['qso:time'] >= instances[instancename]['timelimit']:
            instances[instancename]['qso:stat'] = 'timeout'
        elif outline.startswith('LP Value'):
            instances[instancename]['qso:lpval'] = outline.split()[2].rstrip(',')
        elif outline.find('Trying mpf with') >= 0:
            instances[instancename]['qso:prec'] = max( int(outline.split()[3]), instances[instancename]['qso:prec'] )
        continue

    if outline.startswith('SoPlex version'):
        instances[instancename]['githash'] = outline.split()[-1].rstrip(']')[0:9]
        if not printedIdentifier:
            printedIdentifier = True
            print('\n')
            print(outline)
        continue

    if outline.startswith('Primal solution infeasible') or outline.startswith('Dual solution infeasible'):
        instances[instancename]['status'] = 'fail (solution infeasible)'
        failflag = True
        continue

    if outline.startswith('real:timelimit'):
        instances[instancename]['timelimit'] = fracttofloat(outline.split()[-1])
        continue

    if outline.startswith('Reduced LP has'):
        instances[instancename]['simpcols'] = int(outline.split()[5])
        instances[instancename]['simprows'] = int(outline.split()[3])
        continue

    if outline.startswith('Statistics'):
        stats = True
        continue

    if not stats:
        continue

    inst = instances[instancename]

    if outline.startswith('SoPlex status'):
        beg = outline.find('[')
        reportstat = outline[beg + 1:outline.find(']')] if beg >= 0 else ''
        if reportstat.startswith('optimal'):
            reportstat = 'optimal'
        elif reportstat == 'time limit reached':
            reportstat = 'timeout'
        elif reportstat == 'iteration limit reached':
            reportstat = 'iterlimit'
        elif reportstat == 'objective limit reached':
            reportstat = 'objlimit'
        elif reportstat not in ('infeasible', 'unbounded'):
            reportstat = 'fail (abort)'
        if not failflag or reportstat != 'optimal':
            inst['status'] = reportstat

    elif outline.startswith('  Objective value'):
        inst['value'] = tofloat(fieldval(outline))

    elif outline.startswith('Objective sense'):
        objsenses[instancename] = -1 if 'maximize' in outline.lower() else 1

    elif outline.startswith('  Columns'):
        inorigcol = True
        if 'cols' not in inst:
            inst['cols'] = int(fieldval(outline))
            if 'simpcols' not in inst:
                inst['simpcols'] = inst['cols']
    elif outline.startswith('  Rows'):
        inorigcol = False
        if 'rows' not in inst:
            inst['rows'] = int(fieldval(outline))
            if 'simprows' not in inst:
                inst['simprows'] = inst['rows']

    elif re.match(r' *boxed :', outline):
        inst['boxedcols'] = int(fieldval(outline))
    elif re.match(r' *lower bound :', outline):
        inst['lbcols'] = int(fieldval(outline))
    elif re.match(r' *upper bound :', outline):
        inst['ubcols'] = int(fieldval(outline))
    elif re.match(r' *free :', outline):
        if inorigcol:
            inst['freecols'] = int(fieldval(outline))
        else:
            inst['freerows'] = int(fieldval(outline))
    elif re.match(r' *equal :', outline):
        inst['equalrows'] = int(fieldval(outline))
    elif re.match(r' *ranged :', outline):
        inst['rangedrows'] = int(fieldval(outline))
    elif re.match(r' *lhs :', outline):
        inst['lhsrows'] = int(fieldval(outline))
    elif re.match(r' *rhs :', outline):
        inst['rhsrows'] = int(fieldval(outline))
    elif outline.startswith('  Nonzeros'):
        inst['nonzeros'] = int(fieldval(outline))
    elif re.match(r' *per column :', outline):
        inst['colnonzeros'] = fracttofloat(fieldval(outline))
    elif re.match(r' *per row :', outline):
        inst['rownonzeros'] = fracttofloat(fieldval(outline))
    elif re.match(r' *sparsity :', outline):
        inst['sparsity'] = fracttofloat(fieldval(outline))
    elif re.match(r' *min\. abs\. value :', outline):
        inst['minabsval'] = fracttofloat(fieldval(outline))
    elif re.match(r' *max\. abs\. value :', outline):
        inst['maxabsval'] = fracttofloat(fieldval(outline))

    elif outline.startswith('  Max/sum bound') or outline.startswith('  Max/sum row'):
        inst['pviol'] = vmax(inst.get('pviol', '-'), tofloat(fieldval(outline)))
    elif outline.startswith('  Max/sum redcost') or outline.startswith('  Max/sum dual'):
        inst['dviol'] = vmax(inst.get('dviol', '-'), tofloat(fieldval(outline)))

    elif outline.startswith('Total time'):
        inst['time'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Reading'):
        inst['readtime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Solving'):
        inst['solvetime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Preprocessing'):
        inst['preproctime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Simplex'):
        inst['simplextime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Synchronization'):
        inst['synctime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Transformation'):
        inst['transformtime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Rational'):
        inst['rationaltime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('  Other'):
        inst['othertime'] = fracttofloat(fieldval(outline))

    elif outline.startswith('Refinements'):
        inst['refs'] = int(fieldval(outline))
    elif outline.startswith('  Stalling'):
        if 'stallrefs' not in inst:
            inst['stallrefs'] = int(fieldval(outline))
    elif outline.startswith('  Pivoting'):
        if 'pivrefs' not in inst:
            inst['pivrefs'] = int(fieldval(outline))
    elif outline.startswith('  Feasibility'):
        if 'feasrefs' not in inst:
            inst['feasrefs'] = int(fieldval(outline))
    elif outline.startswith('  Unboundedness'):
        if 'unbdrefs' not in inst:
            inst['unbdrefs'] = int(fieldval(outline))

    elif re.match(r'Iterations +:', outline):
        inst['iters'] = int(fieldval(outline))
        inst['speed'] = round(fracttofloat(inst['iters']) / max(inst.get('solvetime', 0.0), tolerance), 2)
    elif outline.startswith('  From scratch'):
        if 'scratchiters' not in inst:
            inst['scratchiters'] = int(fieldval(outline))
    elif outline.startswith('  From basis'):
        if 'basisiters' not in inst:
            inst['basisiters'] = int(fieldval(outline))
    elif re.match(r'  Primal +:', outline):
        if 'primaliters' not in inst:
            inst['primaliters'] = int(fieldval(outline))
    elif re.match(r'  Dual +:', outline):
        if 'dualiters' not in inst:
            inst['dualiters'] = int(fieldval(outline))
    elif outline.startswith('  Bound flips'):
        if 'flips' not in inst:
            inst['flips'] = int(fieldval(outline))
    elif outline.startswith('  Sol. polishing'):
        if 'polish' not in inst:
            inst['polish'] = int(fieldval(outline))

    elif outline.startswith('LU factorizations'):
        inst['lufacts'] = int(fieldval(outline))
    elif outline.startswith('  Factor. time'):
        inst['factortime'] = fracttofloat(fieldval(outline))
    elif outline.startswith('LU solves'):
        inst['lusolves'] = int(fieldval(outline))
    elif outline.startswith('  Solve time'):
        inst['lusolvetime'] = fracttofloat(fieldval(outline))

    elif decomp:
        if outline.startswith('  Primal Pivots'):
            inst['ppiv'] = int(fieldval(outline))
        elif outline.startswith('  Dual Pivots'):
            inst['dpiv'] = int(fieldval(outline))
        elif outline.startswith('  Primal Candidates'):
            inst['primcand'] = int(fieldval(outline))
        elif outline.startswith('  Dual Candidates'):
            inst['dualcand'] = int(fieldval(outline))
        elif outline.startswith('  Average Primal'):
            inst['avgpdegen'] = fracttofloat(fieldval(outline))
        elif outline.startswith('  Average Dual'):
            inst['avgddegen'] = fracttofloat(fieldval(outline))
        elif outline.startswith('Algorithm Iterations'):
            inst['algiter'] = int(fieldval(outline))
        elif outline.startswith('  Reduced Problem'):
            inst['rediter'] = int(fieldval(outline))
        elif outline.startswith('  Comp. Problem'):
            inst['compiter'] = int(fieldval(outline))
        elif outline.startswith('Red. Problem Status'):
            inst['redprob'] = int(fieldval(outline))
        elif outline.startswith('Comp. Problem Status'):
            inst['compprob'] = int(fieldval(outline))

for inst in instances.values():
    for key in ['pviol','dviol']:
        if isinstance(inst.get(key), float):
            inst[key] = float("{:.2e}".format(inst[key]))

# try parsing solution file
check_solu = bool(soluname)

if check_solu:
    solufile = open(soluname,'r')
    for soluline in solufile:
        solu = soluline.split()
        tag = solu[0]
        name = solu[1]
        if tag == '=inf=':
            value = 'infeasible'
        elif tag == '=unbd=':
            value = 'unbounded'
        elif len(solu) >= 3 and typeofvalue(solu[2]) in [int, float]:
            value = fracttofloat(solu[2])
        else:
            value = 'unknown'
        if name in instances:
            instances[name]['soluval'] = value
            status = instances[name]['status']
            if not status.startswith('fail') and status not in limitstats:
                # check solution status
                if value in ['infeasible', 'unbounded']:
                    instances[name]['status'] = 'ok' if status == value else 'fail (solution on infeasible instance)'
                elif value == 'unknown':
                    instances[name]['status'] = 'solved not verified'
                elif status == 'optimal':
                    solverval = instances[name]['value']
                    valtol = reltol * max(abs(solverval), 1.0)
                    diff = objsenses.get(name, 1) * (solverval - value)
                    instances[name]['status'] = 'fail (objective value)' if tag in ('=opt=', '=best=') and diff > valtol or tag in ('=opt=', '=bestdual=') and diff < -valtol else 'ok'
                else:
                    instances[name]['status'] = 'fail (objective value)'
    solufile.close()

# declare instances without reference unverified
for name in instances:
    if instances[name]['status'] in ('optimal', 'infeasible', 'unbounded'):
        instances[name]['status'] = 'solved not verified'

# save dictionary to file later use in compare script
with open(dataname, 'w') as f:
    json.dump(instances, f)

# count solution status while determining total and (shifted) geometric mean of iterations and time
passed = 0
limits = 0
fails = 0
itergeomshift = 100.0
timegeomshift = 1.0
stotiters = 0
stottime = 0.0
itergeom = 0.0
timegeom = 0.0
shifteditergeom = itergeomshift
shiftedtimegeom = timegeomshift
nprobs = 0
for name in instances:
    status = instances[name]['status']
    if status == 'ok' or status == 'solved not verified':
        passed += 1
    elif status in limitstats:
        limits += 1
    elif status.startswith('fail'):
        fails += 1
    n = instances[name].get('iters') or 0
    t = instances[name].get('time') or 0.0
    stotiters += n
    stottime += t
    nprobs += 1
    itergeom = itergeom ** ((nprobs - 1) / nprobs) * max(n, 1.0) ** (1.0 / nprobs)
    timegeom = timegeom ** ((nprobs - 1) / nprobs) * max(t, 1.0) ** (1.0 / nprobs)
    shifteditergeom = shifteditergeom ** ((nprobs - 1) / nprobs) * max(n + itergeomshift, 1.0) ** (1.0 / nprobs)
    shiftedtimegeom = shiftedtimegeom ** ((nprobs - 1) / nprobs) * max(t + timegeomshift, 1.0) ** (1.0 / nprobs)
shifteditergeom -= itergeomshift
shiftedtimegeom -= timegeomshift

# specify columns for the output (can be modified)
columns = ['name','rows','cols','simprows','simpcols','pviol','dviol','iters','polish','refs','solvetime','value']
if printperf:
    columns.extend(perfcolumns)
columns.append('status')
ncolumns = len(columns)
if decomp:
    columns.extend(decompcolumns)

length = []
output = ''
# calculate maximum width of each column
for i,c in enumerate(columns):
    length.append(len(c))
    for name in instances:
        length[i] = max(length[i],len(str(instances[name].get(c,''))))
    if i == ncolumns:
        output = output + '  |'
    output = output + ' ' + c.rjust(length[i] + 1)

# print column header
print(output)
print('-'*len(output))

# print data for all instances with the computed length
for name in sorted(instances):
    output = ''
    for i,c in enumerate(columns):
        if i == ncolumns:
            output = output + '  |'
        output = output + ' ' + str(instances[name].get(c, '--')).rjust(length[i] + 1)
    print(output)

print('\nResults (testset '+testset+', settings '+settings+'):')
print('{} total: {} pass, {} limit, {} fail'.format(nprobs, passed, limits, fails))
print('iters total {}, geom. mean {:.3f}, shifted geom. mean {:.3f}'.format(stotiters, itergeom, shifteditergeom))
print('time total {:.3f}, geom. mean {:.3f}, shifted geom. mean {:.3f}'.format(stottime, timegeom, shiftedtimegeom))

# try to check for missing files
check_test = bool(testname)

if not check_test:
    print('No testset file found to check run for completeness.')

if not check_solu:
    print('No solution file found to check objective values.')

if check_test:
    testfile = open(testname,'r')
    printedMissing = False
    for testline in testfile:
        linesplit = testline.split('/')
        linesplit = linesplit[len(linesplit) - 1].rstrip(' \n').rstrip('.gz').rstrip('.GZ').rstrip('.z').rstrip('.Z')
        linesplit = linesplit.split('.')
        instancename = linesplit[0]
        for i in range(1, len(linesplit)-1):
            instancename = instancename + '.' + linesplit[i]
        if not instancename in instances:
            if not printedMissing:
                print('\n')
            print('missing instance: '+instancename)
            printedMissing = True

    testfile.close()
