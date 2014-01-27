#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
import math
import json

# This script compares several SoPlex .json files.
# The first argument is the default run that is compared to all other runs. 

# compute speed-up or slow-down factors for all instances of two settings
# compareValue is the name of the value of the loaded dictionary
# compareName is the new name to be used as storage for the factors
def computeFactor(compareValue, defaultSetting, compareSetting):
    "compute speed-up or slow-down factors for all instances of two settings"
    for instance in results[defaultSetting]:
        value = results[compareSetting][instance][compareValue]
        defvalue = results[defaultSetting][instance][compareValue]
        if defaultSetting == compareSetting:
            factors[compareValue][compareSetting][instance] = defvalue
        elif defvalue == 0:
            factors[compareValue][compareSetting][instance] = 0.0
        else:
            factors[compareValue][compareSetting][instance] = float(value) / float(defvalue)

runs = len(sys.argv)
if runs < 2:
    print 'usage: '+sys.argv[0]+' <soplex_test_run1>.json [<soplex_test_run2>.json ...]'
    quit()

# set compare values, used to identify the values from computeFactor
compareValues = ['time','iters']

# parse testset from first file
testset = sys.argv[1].split('/')[-1].split('.')[1]

results = {}
factors = {}
settings = []

# parse given .json files
for run in range(1,runs):
    dataname = sys.argv[run]
    setting = sys.argv[run].split('/')[-1].split('.')[-2]
    if setting in settings:
        setting = setting+str(run)
    settings.append(setting)

    # check for identical testset
    if not testset == dataname.split('/')[-1].split('.')[1]:
        print 'inconsistent testsets'
        quit()

    with open(dataname) as f:
        results[setting] = json.load(f)

# set default setting to compare with
default = settings[0]

# compute all the comparison factors
for c in compareValues:
    factors[c] = {}
    for s in settings:
        factors[c][s] = {}
        computeFactor(c, default, s)

# extract instance names
instances = results[default].keys()
namelength = 5
for i in instances:
    namelength = max(len(i), namelength)
    
namelength = namelength + 1
timelength = 4
iterlength = 5
for s in settings:
    for i in instances:
        timelength = max(len(str(results[s][i]['time'])), timelength)
        iterlength = max(len(str(results[s][i]['iters'])), iterlength)

timelength = timelength + 2
iterlength = iterlength + 2
factorlength = 6

output1 = ' '*namelength + default.center(iterlength + timelength)
output2 = 'name'.ljust(namelength) + 'iters'.rjust(iterlength) + 'time'.rjust(timelength)
for s in settings[1:]:
    output1 = output1 + ' |' + s.center(iterlength + timelength + 2*factorlength)
    output2 = output2 + ' |' + 'iters'.rjust(iterlength) + 'time'.rjust(timelength)
    output2 = output2 + 'iterQ'.rjust(factorlength) + 'timeQ'.rjust(factorlength)

output = '-'*(namelength) + '+' + '-'*(iterlength + timelength) + '+'
output = output + '-'*(iterlength + timelength + 2*factorlength + 1) + '+'
output = output + '-'*(iterlength + timelength + 2*factorlength)
print output
print output1
print output2
print output

sumtime = [0] * len(settings)
sumiter = [0] * len(settings)

# print data for all instances with the computed length
for i in sorted(instances):
    sumtime[0] = sumtime[0] + results[default][i]['time']
    sumiter[0] = sumiter[0] + results[default][i]['iters']
    output = i.ljust(namelength)
    # print results of default settings
    output = output + '{0:{width}d}'.format(results[default][i]['iters'], width=iterlength)
    output = output + '{0:{width}.2f}'.format(results[default][i]['time'], width=timelength)
    # print results of remaining settings
    for idx, s in enumerate(settings[1:]):
        sumiter[idx+1] = sumiter[idx+1] + results[s][i]['iters']
        sumtime[idx+1] = sumtime[idx+1] + results[s][i]['time']
        output = output + '{0:{width}d}'.format(results[s][i]['iters'], width=iterlength+2)
        output = output + '{0:{width}.2f}'.format(results[s][i]['time'], width=timelength)
        output = output + '{0:{width}.2f}'.format(factors['iters'][s][i], width=factorlength)
        output = output + '{0:{width}.2f}'.format(factors['time'][s][i], width=factorlength)
    print output

print
output = 'sum:'.ljust(namelength) + '{0:{width}d}'.format(sumiter[0], width=iterlength)
output = output + '{0:{width}.2f}'.format(sumtime[0], width=timelength)
for idx, s in enumerate(settings[1:]):
    output = output + '{0:{width}d}'.format(sumiter[idx+1], width=iterlength + 2)
    output = output + '{0:{width}.2f}'.format(sumtime[idx+1], width=timelength) + ' '*2*factorlength
print output

