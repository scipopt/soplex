#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
import math
import json

# todo: incorporate timeouts and inconsistent solution status in general

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

# update the current (shifted) geometric metric
def updateGeoMean(new, mean, count, shift):
    assert mean > 0
    if shift == 0:
        shift = 0.0000001
    return math.pow(float(mean+shift), float(count-1)/float(count)) * math.pow(float(new)+shift, 1.0/float(count)) - shift

# print header lines
def printHeader():
    border = '-'*(namelength) + '+' + '-'*(iterlength + timelength) + '+'
    output1 = '-'.join([version[0],hash[0],default]).rjust(namelength + iterlength + timelength)
    output2 = 'name'.ljust(namelength) + '|' + 'iters'.rjust(iterlength-1) + 'time'.rjust(timelength)
    for i,s in enumerate(settings[1:]):
        output1 = output1 + ' |' + '-'.join([version[i+1],hash[i+1],s]).center(iterlength + timelength + 2*factorlength)
        output2 = output2 + ' |' + 'iters'.rjust(iterlength) + 'time'.rjust(timelength)
        output2 = output2 + 'iterQ'.rjust(factorlength) + 'timeQ'.rjust(factorlength)
        border = border + '-'*(iterlength + timelength + 2*factorlength + 1) + '+'
    print border
    print output1
    print output2
    print border

if len(sys.argv) < 2:
    print 'compare several runs of the same testset'
    print 'usage: '+sys.argv[0]+' [ignore=<instance1>,<instance2>,...] <soplex_test_run1>.json [<soplex_test_run2>.json ...]'
    quit()

# set compare values, used to identify the values from computeFactor
compareValues = ['solvetime','iters']

# look for instances to ignore
ignore = []
for a in sys.argv[1:]:
    if a.startswith('ignore='):
        ignore = a.lstrip('ignore=').split(',')
        sys.argv.remove(a)

# parse testset from first file
testset = sys.argv[1].split('/')[-1].split('.')[1]

runs = len(sys.argv)
results = {}
factors = {}
settings = []
version = []
opt = []
hash = []

# load given .json files
for run in range(1,runs):
    dataname = sys.argv[run]
    setting = sys.argv[run].split('/')[-1].split('.')[-2]
    if setting in settings:
        setting = setting+str(run)
    settings.append(setting)
    version.append('.'.join(sys.argv[run].split('/')[-1].split('.')[2:-6]).lstrip('soplex-'))
    opt.append(sys.argv[run].split('/')[-1].split('.')[-3])
    # check for identical testset
    if not testset == dataname.split('/')[-1].split('.')[1]:
        print 'inconsistent testsets'
        quit()
    with open(dataname) as f:
        results[setting] = json.load(f)

# set default setting to compare with
default = settings[0]

# extract instance names
instances = results[default].keys()
namelength = 12
for i in instances:
    namelength = max(len(i), namelength)

# extract git hashes
for s in settings:
    hash.append(results[s][instances[0]]['hash'])

# check all settings for aborts or instances to ignore and remove them
aborts = ''
for s in settings:
    for i in instances:
        if results[s][i]['status'] == 'abort' or i in ignore:
            aborts = aborts + i + '\n'
            instances.remove(i)
            del results[s][i]

# compute all the comparison factors
for c in compareValues:
    factors[c] = {}
    for s in settings:
        factors[c][s] = {}
        computeFactor(c, default, s)

namelength = namelength + 1
timelength = 6
iterlength = 6
for s in settings:
    for i in instances:
        timelength = max(len(str(results[s][i]['solvetime'])), timelength)
        iterlength = max(len(str(results[s][i]['iters'])), iterlength)

timelength = timelength + 2
iterlength = iterlength + 2
factorlength = 8

printHeader()

sumtime = [0] * len(settings)
sumiter = [0] * len(settings)
meantime = [1.0] * len(settings)
meaniter = [1.0] * len(settings)
shmeantime = [1.0] * len(settings)
shmeaniter = [1.0] * len(settings)
shiftiter = 10
shifttime = 0.1
count = 0

# print data for all instances with the computed length
for i in sorted(instances):
    count = count + 1
    time = results[default][i]['solvetime']
    iter = results[default][i]['iters']
    sumtime[0] = sumtime[0] + time
    sumiter[0] = sumiter[0] + iter
    meantime[0] = updateGeoMean(time, meantime[0], count, 0)
    meaniter[0] = updateGeoMean(iter, meaniter[0], count, 0)
    shmeantime[0] = updateGeoMean(time, shmeantime[0], count, shifttime)
    shmeaniter[0] = updateGeoMean(iter, shmeaniter[0], count, shiftiter)
    output = i.ljust(namelength)
    # print results of default settings
    output = output + '{0:{width}d}'.format(iter, width=iterlength)
    output = output + '{0:{width}.2f}'.format(time, width=timelength)
    # print results of remaining settings
    for idx, s in enumerate(settings[1:]):
        time = results[s][i]['solvetime']
        iter = results[s][i]['iters']
        sumiter[idx+1] = sumiter[idx+1] + iter
        sumtime[idx+1] = sumtime[idx+1] + time
        meantime[idx+1] = updateGeoMean(time, meantime[idx+1], count, 0)
        meaniter[idx+1] = updateGeoMean(iter, meaniter[idx+1], count, 0)
        shmeantime[idx+1] = updateGeoMean(time, shmeantime[idx+1], count, shifttime)
        shmeaniter[idx+1] = updateGeoMean(iter, shmeaniter[idx+1], count, shiftiter)
        output = output + '{0:{width}d}'.format(iter, width=iterlength+2)
        output = output + '{0:{width}.2f}'.format(time, width=timelength)
        output = output + '{0:{width}.2f}'.format(factors['iters'][s][i], width=factorlength)
        output = output + '{0:{width}.2f}'.format(factors['solvetime'][s][i], width=factorlength)
    print output

printHeader()

# print summary of comparision
output1 = 'sum:'.ljust(namelength) + '{0:{width}d}'.format(sumiter[0], width=iterlength)
output1 = output1 + '{0:{width}.2f}'.format(sumtime[0], width=timelength)
output2 = 'geo mean:'.ljust(namelength) + '{0:{width}.2f}'.format(meaniter[0], width=iterlength)
output2 = output2 + '{0:{width}.2f}'.format(meantime[0], width=timelength)
output3 = 'shifted:'.ljust(namelength) + '{0:{width}.2f}'.format(shmeaniter[0], width=iterlength)
output3 = output3 + '{0:{width}.2f}'.format(shmeantime[0], width=timelength)
for idx, s in enumerate(settings[1:]):
    output1 = output1 + '{0:{width}d}'.format(sumiter[idx+1], width=iterlength + 2)
    output1 = output1 + '{0:{width}.2f}'.format(sumtime[idx+1], width=timelength) + ' '*2*factorlength
    output2 = output2 + '{0:{width}.2f}'.format(meaniter[idx+1], width=iterlength + 2)
    output2 = output2 + '{0:{width}.2f}'.format(meantime[idx+1], width=timelength) + ' '*2*factorlength
    output3 = output3 + '{0:{width}.2f}'.format(shmeaniter[idx+1], width=iterlength + 2)
    output3 = output3 + '{0:{width}.2f}'.format(shmeantime[idx+1], width=timelength) + ' '*2*factorlength
print output1
print output2
print output3

# print aborted and ignored instances
if not aborts == '':
    print '\naborted and ignored instances:'
    print aborts
