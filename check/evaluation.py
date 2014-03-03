#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
import math
import json

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

    return str

if not len(sys.argv) == 2:
    print 'usage: '+sys.argv[0]+' <soplex_test_run>.out'
    quit()

# specify columns for the output (can be modified)
columns = ['rows','cols','primalviol','dualviol','iters','flips','solvetime','value','status']

outname = sys.argv[1]
dataname = outname.replace('.out','.json')

outfile = open(outname,'r')
outlines = outfile.readlines()
outfile.close()

testset = outname.split('/')[-1].split('.')[1]
testname = 'testset/'+testset+'.test'
soluname = 'testset/'+testset+'.solu'

# maximum length of instance names
namelength = 18

# tolerance for solution value check
tolerance = 1e-6

instances = {}
stats = False
printedIdentifier = False

for idx, outline in enumerate(outlines):
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
            instancename = instancename[length-namelength-2:length-2]

        # initialize new data set
        instances[instancename] = {}
        instances[instancename]['status'] = 'abort'
        # wait for statistics block
        stats = False
    elif outline.startswith('SoPlex version'):
        instances[instancename]['hash'] = outline.split()[-1].rstrip(']')[0:9]
        if not printedIdentifier:
            printedIdentifier = True
            print
            print outline

    # invalidate instancename
    elif outline.startswith('=ready='):
        instancename = ''

    elif outline.startswith('Primal solution infeasible') or outline.startswith('Dual solution infeasible'):
        instances[instancename]['status'] = 'fail'

    elif outline.startswith('Statistics'):
        stats = True

    if stats:
        if outline.startswith('SoPlex status') and not instances[instancename]['status'] == 'fail':
            if outline.find('time limit') >= 0:
                instances[instancename]['status'] = 'timeout'
            else:
                instances[instancename]['status'] = outline.split()[-1].strip('[]')

        elif outline.startswith('Solution'):
            instances[instancename]['value'] = float(outlines[idx+1].split()[-1])

        elif outline.startswith('Original problem'):
            instances[instancename]['cols'] = int(outlines[idx+1].split()[2])
            instances[instancename]['boxedcols'] = int(outlines[idx+2].split()[2])
            instances[instancename]['lbcols'] = int(outlines[idx+3].split()[3])
            instances[instancename]['ubcols'] = int(outlines[idx+4].split()[3])
            instances[instancename]['freecols'] = int(outlines[idx+5].split()[2])
            instances[instancename]['rows'] = int(outlines[idx+6].split()[2])
            instances[instancename]['rangedrows'] = int(outlines[idx+7].split()[2])
            instances[instancename]['lhsrows'] = int(outlines[idx+8].split()[2])
            instances[instancename]['rhsrows'] = int(outlines[idx+9].split()[2])
            instances[instancename]['freerows'] = int(outlines[idx+10].split()[2])
            instances[instancename]['nonzeros'] = int(outlines[idx+11].split()[2])
            instances[instancename]['colnonzeros'] = float(outlines[idx+12].split()[3])
            instances[instancename]['rownonzeros'] = float(outlines[idx+13].split()[3])
            instances[instancename]['sparsity'] = float(outlines[idx+14].split()[2])
            instances[instancename]['minabsval'] = float(outlines[idx+15].split()[4])
            instances[instancename]['maxabsval'] = float(outlines[idx+16].split()[4])

        elif outline.startswith('Iterations'):
            instances[instancename]['iters'] = int(outline.split()[2])
            instances[instancename]['scratchiters'] = int(outlines[idx+1].split()[3])
            instances[instancename]['basisiters'] = int(outlines[idx+2].split()[3])
            instances[instancename]['primaliters'] = int(outlines[idx+3].split()[2])
            instances[instancename]['dualiters'] = int(outlines[idx+4].split()[2])
            instances[instancename]['flips'] = int(outlines[idx+5].split()[3])

        elif outline.startswith('Total time'):
            instances[instancename]['time'] = float(outline.split()[3])
            instances[instancename]['readtime'] = float(outlines[idx+1].split()[2])
            instances[instancename]['solvetime'] = float(outlines[idx+2].split()[2])
            instances[instancename]['preproctime'] = float(outlines[idx+3].split()[2])
            instances[instancename]['simplextime'] = float(outlines[idx+4].split()[2])
            instances[instancename]['synctime'] = float(outlines[idx+5].split()[2])
            instances[instancename]['transformtime'] = float(outlines[idx+6].split()[2])
            instances[instancename]['othertime'] = float(outlines[idx+7].split()[2])

        elif outline.startswith('Refinements'):
            instances[instancename]['refinements'] = int(outline.split()[2])
            instances[instancename]['stalling'] = int(outlines[idx+1].split()[2])

        elif outline.startswith('Violation'):
            primviol = outlines[idx+2].split()[3]
            dualviol = outlines[idx+3].split()[3]
            if typeofvalue(primviol) in [int,float] and typeofvalue(dualviol) in [int,float]:
                instances[instancename]['primalviol'] = float(primviol)
                instances[instancename]['dualviol'] = float(dualviol)
            else:
                instances[instancename]['primalviol'] = '-'
                instances[instancename]['dualviol'] = '-'

        elif outline.startswith('LU factorizations'):
            instances[instancename]['lufacts'] = int(outline.split()[3])
            instances[instancename]['factortime'] = float(outlines[idx+2].split()[3])

        elif outline.startswith('LU solves'):
            instances[instancename]['lusolves'] = int(outline.split()[3])
            instances[instancename]['lusolvetime'] = float(outlines[idx+2].split()[3])

# try parsing solution file
check_solu = False
try:
    with open(soluname):
        check_solu = True
except IOError:
    check_solu = False

if check_solu:
    solufile = open(soluname,'r')
    for soluline in solufile:
        solu = soluline.split()
        tag = solu[0]
        name = solu[1]
        if len(solu) == 3:
            value = solu[2]
            if typeofvalue(value) in [int,float]:
                value = float(value)
        else:
            if tag == '=inf=':
                value = 'infeasible'
            else:
                value = 'unknown'
        if name in instances:
            instances[name]['soluval'] = value
            if not instances[name]['status'] in ['timeout', 'fail', 'abort']:
                # check solution status
                if value in ['infeasible', 'unbounded']:
                    if not instances[name]['status'] == value:
                        instances[name]['status'] = 'fail'
                elif instances[name]['value'] == 'optimal':
                    if (abs(instances[name]['value'] - value))/max(abs(instances[name]['value']),abs(value)) > tolerance:
                        instances[name]['status'] = 'inconsistent'
    solufile.close()

# save dictionary to file later use in compare script
with open(dataname, 'w') as f:
    json.dump(instances, f)

# count solution status
fails = sum(1 for name in instances if instances[name]['status'] == 'fail')
timeouts = sum(1 for name in instances if instances[name]['status'] == 'timeout')
infeasible = sum(1 for name in instances if instances[name]['status'] == 'infeasible')
optimal = sum(1 for name in instances if instances[name]['status'] == 'optimal')
aborts = sum(1 for name in instances if instances[name]['status'] == 'abort')

length = []

output = 'name'.ljust(namelength)
# calculate maximum width of each column
for i,c in enumerate(columns):
    length.append(len(c))
    for name in instances:
        length[i] = max(length[i],len(str(instances[name].get(c,''))))
    output = output + ' ' + c.rjust(length[i] + 1)

# print column header
print output
print '-'*len(output)

# print data for all instances with the computed length
for name in sorted(instances):
    output = name.ljust(namelength)
    for i,c in enumerate(columns):
        output = output + ' ' + str(instances[name].get(c, '--')).rjust(length[i] + 1)
    print output

print
print 'Results (testset '+testname.split('/')[-1].split('.')[-2]+', settings '+outname.split('/')[-1].split('.')[-2]+'):'
print '{} total: {} optimal, {} fails, {} timeouts, {} infeasible, {} aborts'.format(len(instances),optimal,fails,timeouts,infeasible,aborts)

# try to check for missing files
check_test = False
try:
    with open(testname):
        check_test = True
except IOError:
    print 'No testset file found to check run for completeness.'

if not check_solu:
    print 'No solution file found to check objective values.'

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
        length = len(instancename)
        if length > namelength:
            instancename = instancename[length-namelength-2:length-2]
        if not instancename in instances:
            if not printedMissing:
                print
            print 'missing instance: '+instancename
            printedMissing = True

    testfile.close()
