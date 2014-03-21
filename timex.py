'''measure runtime of script - just start and stops clock and calculate stop-start'''
__author__ = 'jmadzo'

import time		#to time how long the script is running
def startTime():
    '''start timer'''
    return time.time(), time.localtime()
    
def stopTime(start,startLocal):
    '''stop timer and print running times'''
    elapsedLocal=time.localtime()
    elapsed=(time.time()-start)
    #print info about runtime
    print;print 'All done, took: ~ ',
    print '-> (',elapsed,'sec.)'
    print
    print 'started  ->',':'.join([str(i) for i in startLocal[3:6]])
    print 'finished ->',':'.join([str(i) for i in elapsedLocal[3:6]])
    return None

def today():
    '''find today's date and return it in format MM-DD-YYYY'''
    import time
    return '-'.join([str(x) for x in time.localtime()[0:3]])