# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:27:05 2016

@author: ab9621
"""

import subprocess
#import shutil

#def SubmitRun(jobname, path):
def SubmitRun(jobname):
    #jobname should be e.g. 'inputfile' NOT 'inputfile.inp'
    #shutil.copyfile('{}.inp'.format(jobname), '{}/{}.inp'.format(path, jobname))
    #subprocess.call('pogoFromAb {}.inp'.format(jobname), cwd=path)
    #subprocess.call('pogoBlock {}.pogo-inp'.format(jobname), cwd=path)
    #subprocess.call('pogoSolve {} -o'.format(jobname), cwd=path)
    subprocess.call('pogoFromAb {}.inp'.format(jobname))
    subprocess.call('pogoBlock {}.pogo-inp'.format(jobname))
    subprocess.call('pogoSolve {} -o'.format(jobname))


if __name__ == '__main__':
    job = 'Attempt1'
    #path = 'C:/Google Drive/pogo'
    SubmitRun(job )