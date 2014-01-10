#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
import os
import math
import pandas as pd
from decimal import Decimal


soluname = sys.argv[1]
testname = sys.argv[2]
outname = sys.argv[3]

solufile = open(soluname,'r')
testfile = open(testname,'r')
outfile = open(outname,'r')

solulines = solufile.readlines()
outlines = outfile.readlines()

# print identifier
print
print outlines[2]

# maximum length of instance names
namelength = 18

# tolerance for solution value check
tolerance = 1e-6

instances = {}
data = {'status': 'unknown',
        'value': '',
        'rows': '',
        'cols': '',
        'presolrows': '',
        'presolcols': '',
        'iters': '',
        'primaliters': '',
        'dualiters': '',
        'time': '',
        'primalviol': '',
        'dualviol': '',
        'solustat': '',
        'soluval': ''}

for idx, outline in enumerate(outlines):
    if outline.startswith('@01'):
        # convert line to used instance name
        linesplit = outline.split('/')
        linesplit = linesplit[len(linesplit) - 1].rstrip(' \n').rstrip('.gz').rstrip('.GZ').rstrip('.z').rstrip('.Z')
        linesplit = linesplit.split('.')
        instancename = linesplit[0]
        for i in range(1, len(linesplit)-1):
            instancename = instancename + '.' + linesplit[i]
        length = len(instancename)
        if length > namelength:
            instancename = instancename[length-namelength-2:length-2]

        # initialize empty entry
        instances[instancename] = {key:val for key,val in data.iteritems()}

    # invalidate instancename
    if outline.startswith('=ready='):
        instancename = ''

    elif outline.startswith('SoPlex status'):
        instances[instancename]['status'] = outline.split()[-1].strip('[]')

    elif outline.startswith('Solution'):
        instances[instancename]['value'] = float(outlines[idx+1].split()[-1])

    elif outline.startswith('Original problem'):
        instances[instancename]['cols'] = int(outlines[idx+1].split()[-1])
        instances[instancename]['rows'] = int(outlines[idx+6].split()[-1])

    elif outline.startswith('Iterations'):
        instances[instancename]['iters'] = int(outline.split()[-1])

    elif outline.startswith('  Primal'):
        instances[instancename]['primaliters'] = int(outline.split()[-2])

    elif outline.startswith('  Dual'):
        instances[instancename]['dualiters'] = int(outline.split()[-2])

    elif outline.startswith('Total time'):
        instances[instancename]['time'] = float(outline.split()[-1])

    elif outline.startswith('Violations'):
        instances[instancename]['primalviol'] = float(outlines[idx+2].split()[-3])
        instances[instancename]['dualviol'] = float(outlines[idx+3].split()[-3])

    elif outline.startswith('Primal solution infeasible') or outline.startswith('Dual solution infeasible'):
        instances[instancename]['status'] = 'fail'
        fail = fail + 1

# parse data from solufile
for soluline in solulines:
    tag, name, value = soluline.split()
    value = float(value)
    if instances.has_key(name):
        instances[name]['soluval'] = value
        # check solution status
        if value in ['infeasible', 'unbounded']:
            if not instances[name]['status'] == value:
                instances[name]['status'] = 'fail'
        else:
            if (abs(instances[name]['value'] - value))/max(abs(instances[name]['value']),abs(value)) > tolerance:
                instances[name]['status'] = 'inconsistent'

df = pd.DataFrame(instances).T

fails = sum(1 for name in instances.keys() if instances[name]['status'] == 'fail')
timeouts = sum(1 for name in instances.keys() if instances[name]['status'] == 'timeout')

print df[['rows','cols','primalviol','dualviol','value','iters','time','status']].to_string(float_format=lambda x:"%6.9g"%(x))
print
print 'total: {} fails: {} timeout: {}'.format(len(instances),fails,timeouts)
