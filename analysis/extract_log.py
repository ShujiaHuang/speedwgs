#!/usr/bin/env python
import os
import json
import math
import datetime
import time
import sys
import signal
import subprocess
import threading
import urllib2

instanceID = sys.argv[1]
taskName = sys.argv[2]
threadNum = int(sys.argv[3])
logToken = sys.argv[4]

logUrl = "http://10.101.207.200:22224/api/projects/hadoop_on_odps/instances/{0}?log&logtype=Stderr&size=10000&id={1}&authorization_token={2}"

def getLogTime(instanceID, logID, reducerID):
    logJsonText = []
    try:
        logJsonText = urllib2.urlopen(logUrl.format(instanceID, logID, logToken)).readlines()
    except urllib2.URLError as e:
        sys.stderr.write(e.reason + '\n')

    commandTime = []
    isBegin = False

    beginHint = ">> reducer begin "
    def checkBegin(line):
        idx = line.find(beginHint)
        return idx != -1 and int(line[idx + len(beginHint):]) == reducerID
    logHint = "<< command finished in "
    endHint = "everything seems good, exit"
    def checkEnd(line):
        return isBegin and line.find(endHint) != -1
    for logLine in logJsonText:
        if not isBegin:
            if not checkBegin(logLine):
                continue
            isBegin = True
        if checkEnd(logLine):
            break
        idx = logLine.find(logHint)
        if idx == -1:
            continue
        commandTime.append(int(logLine[idx + len(logHint):].split(' ', 1)[0]))
    if len(commandTime) == 0:
        return [0] * 8 # Should be 8 numbers if everything is right
    return commandTime

jobJsonCmd = 'odpscmd -e "http GET /projects/hadoop_on_odps/instances/{0}?instancedetail&taskname={1}"'.format(instanceID, taskName)
jobJsonText = os.popen(jobJsonCmd).readlines()

# omit http response header
while jobJsonText[0][0] != '{':
    jobJsonText = jobJsonText[1:]
jobJsonData = json.loads("".join(jobJsonText))
reducers = jobJsonData['Instance']['Stages'][1]['Workers']
reducerStat = {} # { reducerID: [totalTime, inputCount] }
print "reducerID", "inputCount", "totalTime", "rmdupTime", "indexTime", "realignTargetCreateTime", "indelRealignTime", "baseRecTime", "printReadsTime", "UniGenotyTime"

def processReducers(reducers):
    for r in reducers:
        if (r['Status'] == 'Interrupted'):
            continue
        reducerID = int(r['WorkerID'].split("/")[2].split("#")[1].split("_")[0])
        totalTime = long(r['EndTime'])-long(r['StartTime'])
        counter = json.loads(r['GblCounter'])
        inputCount = long(counter['InputCounters']['input']['records'])
        rmdup_t, index_t, javaLn_t, realignTargetCreate_t, indelRealign_t, baseRec_t, printReads_t, uniGenoty_t = getLogTime(instanceID, r['LogID'], reducerID)
        sys.stderr.write(str(reducerID) + '\n')
        sys.stdout.write(' '.join(map(str, [reducerID, inputCount, totalTime, rmdup_t, index_t, realignTargetCreate_t, indelRealign_t, baseRec_t, printReads_t, uniGenoty_t])) + '\n')

redPerThread = int(math.ceil(1.0 * len(reducers) / threadNum))
threads = []
for i in range(0, threadNum):
    try:
        t = threading.Thread(target=processReducers, args=(reducers[i * redPerThread: (i + 1) * redPerThread],))
        threads.append(t)
        t.start()
    except:
        print "ERROR, unable to start"

for t in threads:
    t.join()
