# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:27:05 2016

@author: ab9621
"""

import os
import subprocess
#import shutil
import time

#def SubmitRun(jobname, path):
def SubmitRun(jobname):
    #jobname should be e.g. 'inputfile' NOT 'inputfile.inp'
    #shutil.copyfile('{}.inp'.format(jobname), '{}/{}.inp'.format(path, jobname))
    #subprocess.call('pogoFromAb {}.inp'.format(jobname), cwd=path)
    #subprocess.call('pogoBlock {}.pogo-inp'.format(jobname), cwd=path)
    #subprocess.call('pogoSolve {} -o'.format(jobname), cwd=path)
    time1 = time.time_ns()
    os.system('cmd /c "abaqus job={} inp={}.inp"'.format(jobname, jobname))
    # subprocess.call('abaqus job={}, inp={}.inp'.format(jobname, jobname))
    time2 = time.time_ns()
    print('{} Runtime: {:1.3g}s'.format(jobname, (time2-time1)*10**-9))


if __name__ == '__main__':
    os.chdir('C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v1\\io files')
    for element in range(1):
        job = '16El_FMC_el_{}'.format(element)
        SubmitRun(job)
        time.sleep(5)