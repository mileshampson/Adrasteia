""" Performance testing for Adastreia
    Copyright (C) 2006   Miles Hampson

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
""" 
import math
import random
import string
import time

import Adrasteia
from Adrasteia import *
from Align_score import *
import GenSequence
      
#------------------------------------------------------------------------------#
#Generates random RSDB-60 protein sequences. The length of each sequence is
#taken from a discretized log-normal distribution, between min_bound and max_bound 
#if it is not 0.The amino acids of a sequence are generated by an independent, 
#identically distributed process.  The probabilities for that distribution are 
#selected from a mixture of Dirichlet densities.
#------------------------------------------------------------------------------#     
def generateRandomSequences(runs,min_bound,max_bound):
    GenSequence.generate(runs)
    f=open('GeneratedSequences.txt', 'r')
    sequences = f.read().split()
    random.seed(time.clock())
    for i in xrange(runs):
        p1 = int(math.floor(len(sequences[i])*random.uniform(0,1)))
        p2 = int(math.floor(len(sequences[i])*random.uniform(0,1)))
        while (abs(p2-p1)>=max_bound) or (abs(p2-p1)<min_bound) or p1==p2:
            p1 = int(math.floor(len(sequences[i])*random.uniform(0,1)))
            p2 = int(math.floor(len(sequences[i])*random.uniform(0,1)))
        if(p2>p1):
            sequences[i] = sequences[i][p1:p2]
        else:
            sequences[i] = sequences[i][p2:p1]
    f.close()
    return sequences
    
#-------------------------------------------------------#
#As above but generates sequences of the specified length
#-------------------------------------------------------#          
def generateRandomSequencesss(runs,length):
    GenSequence.generate(runs)
    f=open('GeneratedSequences.txt', 'r')
    sequences = f.read().split()
    random.seed(time.clock())
    for i in xrange(runs):
        p1 = int(math.floor(len(sequences[i])*random.uniform(0,1)))
        p2 = p1+length
        sequences[i] = sequences[i][p1:p2]
    return sequences
    
    
def generateBackgroundProfile(file_to_parse,test_num):
    profile=open('BackgroundProfile.txt', 'w')
    db_file = open(file_to_parse, 'r')    
    means = []
    stdevs = []
    for i in range(1,80):
        seqs = generateRandomSequences(test_num,i)        
        cum_align_score = 0
        cum_square_score = 0
        trials = 0
        for j in range(test_num):
            record = []
            db_file.seek(0)
            for line in db_file:           
                #database body section
                if line[0] != '>':
                    record.append(line[:-1])               
                #database header
                else:
                    record = ''.join(record)
                    segments = [record[k+(i/2):k+i] for k in range(0,len(record),i*2)]
                    for l in range(len(segments)):
                        trials += 1
                        al_score = score(i,seqs[j],segments[l])
                        cum_align_score += al_score
                        cum_square_score += pow(al_score,2)                       
                    record = []
        mean = float(cum_align_score)/float(trials)
        means.append(mean)
        stdevs.append(pow(fabs((float(cum_square_score)/float(trials))-pow(mean,2)),2))
        print "done %i" % (i)
    profile.write("means = [",)
    for i in range(1,80):
        if i==79:
            profile.write("%.4f" % (means[i-1]))
        else:
            profile.write("%.4f, " % (means[i-1]))
    profile.write("]\nstdevs = [",) 
    for i in range(1,80):   
        if i==79:
            profile.write("%.4f" % (stdevs[i-1]))
        else:
            profile.write("%.4f, " % (stdevs[i-1]))
    profile.write("]\n")
    db_file.close()
    profile.close()
            
    

#---------------------------------------------------------------------------------#
#Print the average time, over the given number of runs, taken to find matches that 
#are edit_distance away from a generated 
#sequence of min < size < max, in an online search of the given database
#---------------------------------------------------------------------------------#    
def timeDirectMatch(runs,min_bound,max_bound,file_to_parse,ed):
    elapsedtime= 0
    i = 0
    seqs = generateRandomSequences(runs,min_bound,max_bound)
    print 'TimeOnlineSequence - measuring matching times...'
    for i in xrange(runs):
        if(runs<15):print '->matching sequence: ',seqs[i]
        t1 = time.clock()
        direct_match(file_to_parse,seqs[i],ed)
        t2 = time.clock()
        elapsedtime = elapsedtime + (t2-t1)
    print '->Matching took an average of %0.3fms' % ((elapsedtime*1000.0)/runs)
    
#---------------------------------------------------------------------------------#
#Print the average time, over the given number of runs, taken to find matches that 
#are edit_distance away from a generated sequence of size < max, in an 
#offline search of the given database using the given fingerprint_size
#---------------------------------------------------------------------------------#    
def timePreprocessAndMatch(runs,min_bound,max_bound,file_to_parse,ed):
    elapsedtime= 0
    i = 0
    seqs = generateRandomSequences(runs,min_bound,max_bound)
    print 'TimeOnlineFilteredSequence - preprocessing database...'
    t1 = time.clock()
    m=FastaDatabase(file_to_parse)
    t2 = time.clock()
    print '->Database preprocessing took %0.3fms' % ((t2-t1)*1000.0)
    print 'TimeOnlineFilteredSequence - measuring matching times...'
    for i in xrange(runs):
        if(runs<15):print '->matching sequence: ',seqs[i]
        t1 = time.clock()
        m.match(seqs[i],ed)
        t2 = time.clock()
        elapsedtime += (t2-t1)
    del m
    print '->Matching took an average of %0.3fms' % ((elapsedtime*1000.0)/runs)
    

#---------------------------------#
#As above but uses multiple queries
#---------------------------------#    
def timeMultipleDirectMatch(num,runs,min_bound,max_bound,file_to_parse,ed):
    elapsedtime= 0
    i = 0
    seqs = generateRandomSequences(runs*num,min_bound,max_bound)
    print 'TimeOnlineSequence - measuring matching times...'
    for i in range(0,runs*num,num):
        param = []
        for j in range(num):
            param.append(seqs[i+j])
        t1 = time.clock()
        direct_match(file_to_parse,param,ed)
        t2 = time.clock()
        elapsedtime = elapsedtime + (t2-t1)
    print '->Matching took an average of %0.3fms' % ((elapsedtime*1000.0)/runs)
    
#---------------------------------#
#As above but uses multiple queries
#---------------------------------# 
def timeMultiplePreprocessAndMatch(num,runs,min_bound,max_bound,file_to_parse,ed):
    elapsedtime= 0
    i = 0
    seqs = generateRandomSequences(runs*num,min_bound,max_bound)
    print 'TimeOnlineFilteredSequence - preprocessing database...'
    t1 = time.clock()
    m=FastaDatabase(file_to_parse)
    t2 = time.clock()
    print '->Database preprocessing took %0.3fms' % ((t2-t1)*1000.0)
    print 'TimeOnlineFilteredSequence - measuring matching times...'
    for i in range(0,runs*num,num):
        param = []
        for j in range(num):
            param.append(seqs[i+j])
        #param2 = []
        #param2.append(param)
        t1 = time.clock()
        m.match(param,ed)
        t2 = time.clock()
        elapsedtime += (t2-t1)
    del m
    print '->Matching took an average of %0.3fms' % ((elapsedtime*1000.0)/runs)

