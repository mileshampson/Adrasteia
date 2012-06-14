#!/usr/bin/env python
#file Adrasteia.py
""" Functions for performing approximate matching on short proteins 
    Copyright (C) 2006  Miles Hampson

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
    
    See changelog.txt in top level directory for latest bugs and changes
"""
from marshal import load,dump
from math import floor
from operator import itemgetter
from string import split

from agrepy import *
from Align_score import *
from Bitwise import *

#The maximum p-values to report (scores with p-values above this will
#not be reported back to the user)
PVAL_THRESHOLD = 0.05

#-----Class that holds information about a FASTA database-------------------#
class FastaDatabase:
    """
    This class holds information about a FASTA database which is given to
    it during initialization. It indexes the database using the method of
    superimposed code words, this index is stored in the accompanying
    Bitwise module (implemented in c). Once initialized, the match function
    can be called with a query, and the indices will be used to quickly locate
    any records in the database that have permutations of index data that
    would allow a match. The candidate records (dropset) are then tested using
    a version of the agrep approximate matching function (Apse), which locates
    any that have approximate matches to the query, within some specified
    error tolerance (edit distance). Each match is then given a score, and a
    p-value indicating how likely it is to not to be just a fluke. The results
    are ranked from most to least likely by p-value, and output to the user.
    """
    #Length of the code words
    WORD_WIDTH = 1024
    #Number of bits to set for the encoding of an attribute
    NUM_SET_BITS = 1
    #Each line of the database will be broken into pieces of these sizes 
    #for preprocessing.
    PIECE_SIZES = [2,3]
    #Table of piece sizes that most efficently fit each segment size
    SEGMENTS = [[],[],[2],[3],[2,2],[3,2],[3,3],[3,2,2],[3,3,2],[3,3,3],\
                [3,3,2,2],[3,3,3,2],[3,3,3,3],[3,3,3,2,2],[3,3,3,3,2],\
                [3,3,3,3,3],[3,3,3,3,2,2],[3,3,3,3,3,2],[3,3,3,3,3,3],\
                [3,3,3,3,3,2,2],[3,3,3,3,3,3,2],[3,3,3,3,3,3,3],\
                [3,3,3,3,3,3,2,2],[3,3,3,3,3,3,3,2],[3,3,3,3,3,3,3,3],\
                [3,3,3,3,3,3,3,2,2],[3,3,3,3,3,3,3,3,2],\
                [3,3,3,3,3,3,3,3,3],[3,3,3,3,3,3,3,3,2,2],\
                [3,3,3,3,3,3,3,3,3,2],[3,3,3,3,3,3,3,3,3,3]]
    #Length of a line of data in this database
    LINE_LENGTH = 60
    #For efficency set a max bound on the size of any query protein
    MAX_QUERY_LENGTH = 60

    #---------------------------------------------------#
    #Constructor. Sets the instance to the given database
    #---------------------------------------------------#
    def __init__(self,db_file_name):
        
        #-------------
        #1. Initialize
        #-------------
        #Check they have specified a string for the database name
        if type(db_file_name) is not type("abc"):
            error = "Error: Invalid type for database name"
            print error
            return error
        #in memory index of the position of each record header
        self.header_index = {}
        #store the name of our database
        self.database_name = db_file_name
        #the name of the index for the database
        self.index_name = db_file_name[0:-6]+".idx"
        #the name of the additional slices file
        self.slice_file_name = db_file_name[0:-6]+".slices"
        #initialise the bitwise module
        initialize(self.WORD_WIDTH,self.NUM_SET_BITS)
        
        #---------------------------
        #2. Unmarshal data from disk
        #---------------------------
        #check if there is a pre-existing index for
        #this db, and if so use the data in there
        try:
            file_index = load(open(self.index_name, "rb"))
            #There is an index, so load data from it
            #The data is in the format {record_segments,
            #bit_string_elms,header_index}
            setDBInfo(file_index[0],file_index[1])
            self.header_index = file_index[2]
            loadBitSlices(0,len(self.PIECE_SIZES)-1,len(self.slice_file_name),
                                              self.slice_file_name+str(0))
            print "Sucessfully loaded index for "+db_file_name
            return None
        except IOError:
            print "No pre-existing index, creating one now."
            print "Please be patient as this takes about one"
            print "minute per Mb. This creates .idx and .slices"
            print "files next to the database so that index"
            print "creation can be skipped for future sessions."
        try:
            db_file = open(db_file_name, 'r')
        except IOError:
            error = "Error: Cannot open specified db file"
            raise Exception(error)
        
        #-------------------
        #3. Process Database
        #-------------------
        #put frequently accessed globals into local vars
        header_index = self.header_index
        bit_string_elms=0
        #Number of different pieces we are going to encode
        num_pieces = len(self.PIECE_SIZES)
        for i in xrange(num_pieces):
            piece_size = self.PIECE_SIZES[i]
            #db index pointer
            file_index_pointer=-1
            db_file.seek(0)
            #find the location of the beginning of the line
            position = db_file.tell()
            while 1: 
                #read in each line individually, rather than the read ahead 
                #in 'for line in db_file' - while this method is slower it 
                #allows us to record the position of each record in an 
                #index into the db
                line=db_file.readline()
                if not line:
                    break
                #if we have a database body section, index it, and break it 
                #into pieces which are encoded and stored
                if line[0] != '>': 
                    length = len(line)-piece_size              
                    for j in xrange(length):
                        encodeDBPiece(line[j:j+piece_size])               
                #Otherwise it is a header, so store its info
                else:
                    file_index_pointer +=1
                    header_index[file_index_pointer] = position
                    allocateRecord()
                #before moving on to next line, get its starting position
                position = db_file.tell()
            #end db loop 
            #do a final bit slice
            bitSlice() 
            #Marshal info about this run
            print "    Writing index for run ",i
            writeBitSlices(self.slice_file_name+str(i))
            print "    Done"
            #reset bitwise for next run if there is one
            if i+1<num_pieces:
                bit_string_elms=resetDBInfo()     
        #Finished indexing so close db file
        db_file.close()
        #Marshal index to index file for this database
        #set up Bitwise with correct info to begin matching
        print "    Resetting Database Info"
        record_segments = setDBInfo(0,bit_string_elms)
        print "    Writing Header Index file"
        file_index = [record_segments,bit_string_elms,header_index]
        try:
            index_file = open(self.index_name,"wb")
            dump(file_index,index_file)  
            index_file.close()
        except IOError:
            error = "Error: Cannot save index to file"
            raise Exception(error) 
        #Also load slices into memory
        print "    Loading Slices into Memory"
        loadBitSlices(0,len(self.PIECE_SIZES)-1,len(self.slice_file_name),
                                              self.slice_file_name+str(0))
        print "Sucessfully created index for "+db_file_name

    #----------------------------------------------------------------------#
    #Find matches in database within specified distance of the given pattern
    #----------------------------------------------------------------------# 
    def match(self,patterns,max_errors):
        """ Parameters:
            
            patterns - a string (i.e. "MKFL")
                    - a single list element (i.e. ["MKFL"])
                    - multiple list elements (i.e. ["LIL","MKFL","CLF"])
                    - grouped list elements (i.e.  [["LIL","MKFL"],"CLF"])
                                            
                    The first two cases are equivalent, and will result in
                    a search for all records that contain the string/element
                    
                    Multiple list elements will result in a search for any
                    records that contain ANY of the elements in the list
                    
                    Elements in a list can also be grouped together. If a
                    group contains only a single element it will not have
                    any effect, but specifying more than one element in a
                    group has the effect that only records that contain ALL
                    elements in a group will be searched for (but if there 
                    is more than one group a record needs to contain only one
                    group to match, in other words elements within a group
                    are ANDed together, and the groups are ORed). Also it
                    is possible that false drops will return some subset
                    of the patterns even if they are not explicitly
                    searched for. The main advantage of grouping elements 
                    together like this is dramatically faster searches.
                    
        max-errors - an integer
                    - a list (i.e [2,3,1])
                    
                    A single integer value will be applied to ALL patterns.
                    
                    If a seperate error rate is desired for each pattern then
                    specify a list with 1 integer per pattern

        """    
        #---------------------------
        #1. Check validity of inputs
        #---------------------------
        #put patterns into uniform list format
        #first deal with a single query specified as a string rather than a list
        tmp_patterns = []
        if type(patterns)==type("abc"):
            tmp_patterns.append(patterns)
        else:
            tmp_patterns=patterns
        #now pull out all groups, storing the info of each
        groups = []   
        patterns_to_match = []
        i=0
        for pat in tmp_patterns:
            if type(pat)==type([]):
                group = []
                for elm in pat:
                    patterns_to_match.append(elm)
                    group.append(i)
                    i+=1
                groups.append(group)  
            else:
                patterns_to_match.append(pat)
                i+=1   
        #add an int for patterns without groups
        if len(groups)==0:
            groups.append([-1])
        #Finally get some data about this list 
        try:
            num_patterns = len(patterns_to_match)
        except TypeError:
            error = "Error:  Invalid type of pattern data"
            raise Exception(error)
        #and check at least one query has been specified
        if num_patterns<1:
            error = "Error:  At least one pattern to find must be specified"
            raise Exception(error)   
        #Now move on to checking max_errors    
        if type(max_errors) is type([]):
            i=0
            for elm in max_errors:
                if elm>8:
                    print "An element of max errors is too high, changing it to 8"
                    elm=8
                if elm<1:
                    print "An element of max errors is too low, changing it to 1"
                    elm=1
                i+=1
            if i>num_patterns:
                error = "Error:  More errors were specified than patterns given"
                raise Exception(error)
            if i<num_patterns:
                error = "Error:  Less errors were specified than patterns given"
                raise Exception(error)
        else:
            if type(max_errors) is type(1):
                if max_errors>8:
                    print "The maximum error rate is to high, changing it to 8"
                    max_errors=8
                if max_errors<1:
                    print "The maximum error rate is too low, changing it to 1"
                    max_errors=1
                errors = max_errors
                max_errors = []
                for i in range(num_patterns):
                    max_errors.append(errors)   
            else:
                error = "Error:  Invalid type for number of errors"
                print error
                return error 

        #------------------------
        #2. Encode Query Patterns
        #------------------------
        matchers = []
        query_pieces = []
        piece_sizes,segments = self.PIECE_SIZES,self.SEGMENTS
        i=0
        for pattern in patterns_to_match:
            #create a pattern object for each query
            matchers.append(compile(pattern,len(pattern),max_errors[i]))
            #check the query length is within acceptable bounds
            pattern_length = len(pattern)      
            #if the pattern is too long simply abort, we don't want to deal
            #with overflow issues and there are better tools out there anyway
            if pattern_length>self.MAX_QUERY_LENGTH:
                error = []
                error.append("Error: This program currently cannot deal with ")
                error.append("proteins longer than ")
                error.append(str(self.MAX_QUERY_LENGTH))
                error.append(" aa, sorry")
                error = ''.join(error)
                clearMatchInfo()
                raise Exception(error)
            #break the query into a number of segments equal to the 
            #maximum number of errors for the query plus one
            num_segs = max_errors[i]+1
            segs_size=pattern_length/num_segs
            #if the pattern is too short then abort all 
            #preprocessed matching and give all input to direct
            #match to deal with and return to the user
            if segs_size<piece_sizes[0]:
                print "A query had too many errors for its size"
                print "Using direct_match instead..."
                clearMatchInfo()
                return direct_match(self.database_name,patterns,max_errors)
            #for each segment encode all pieces
            createSegments(num_segs)
            completed=0
            for j in range(num_segs):
                #on the last one make size the leftover
                if j+1==num_segs:
                    segs_size = pattern_length-(segs_size*(num_segs-1))
                for elm in segments[segs_size]:
                    encodeQueryPiece(j,pattern[completed:completed+elm],elm-2)
                    completed+=elm
            for elm in groups:
                if i in elm:
                    updateDropsetAnd()
                else:
                    updateDropsetOr()
            i+=1

        #-------------------------
        #3. Set up database access
        #-------------------------  
        try:
            db_file = open(self.database_name, 'r')
        except IOError:
            error = "Error: Cannot open specified db file"
            raise Exception(error)  
            
        #----------------------------------------------
        #4. Search dropped database records for matches
        #----------------------------------------------   
        #storage for information about each match
        #this allows for later calculation of p values and
        #output of this information to the user
        results_record_headers = []
        results_surrounding_text = []
        num_match_recs = 0
        #handles for called functions
        index = self.header_index
        db_pos = db_file.seek
        #start parsing dropfile
        rec_num=nextDrop(1)
        while(rec_num):  
            db_pos(index[rec_num-1])
            title =  db_file.readline().split()[1]
            rec = []
            line = db_file.readline()
            while(line and line[0]!='>'):
                rec.append(line[:-1])
                line = db_file.readline()
            record = ''.join(rec)
            begin_rec = 0
            for i in range(num_patterns):
                matches = agrepy(patterns_to_match[i], len(patterns_to_match[i]),\
                                 record,len(record),1,matchers[i])
                if matches:
                    if begin_rec==0:
                        begin_rec = 1
                        num_match_recs +=1
                    results_record_headers.append(title)
                    pattern_matches = []
                    for elm in matches:
                        start=0
                        if elm[0]>0:
                            start = elm[0]
                        pattern_matches.append(record[start:elm[1]]) 
                    results_surrounding_text.append(pattern_matches) 
            rec_num=nextDrop(rec_num+1)
        clearMatchInfo()      
        db_file.close()
        
        #--------------------
        #5. Deal with results
        #-------------------- 
        outputResults(patterns_to_match,num_match_recs,results_record_headers,
                      results_surrounding_text)
        #Finally return result headers to the caller in case
        #they wish to make use of them in some other way
        return results_record_headers

    #-----------------------
    # Free memory in Bitwise
    #----------------------- 
    def __del__(self):
        clearDatabase()

#------------------------------------------------------------------------------#
#Find matches in database within specified distance of the given pattern
#This version does not require a preprocessed database, any file in FASTA
#format can be specified as the database parameter.
#------------------------------------------------------------------------------#   
def direct_match(db_file_name,patterns,max_errors):
    try:
        num_patterns = len(patterns)
    except TypeError:
        error = "Error:  Invalid type of pattern data"
        raise Exception(error) 
    if type(max_errors) is type([1]):
        for elm in max_errors:
            if elm>8:
                print "An element of max errors is too high, changing it to 8"
                elm=8
            if elm<1:
                print "An element of max errors is too low, changing it to 1"
                elm=1
    else:
        if type(max_errors) is type(1):
            if max_errors>8:
                print "The maximum error rate is to high, changing it to 8"
                max_errors=8
            if max_errors<1:
                print "The maximum error rate is too low, changing it to 1"
                max_errors=1
            errors = max_errors
            max_errors = []
            for i in range(num_patterns):
                max_errors.append(errors)  
        else:
            error = "Error:  Invalid type for number of errors"
            raise Exception(error) 
    #check at least one query has been specified
    if num_patterns<1:
        error = "Error: At least one pattern to find must be specified"
        raise Exception(error)
    matchers = []
    #deal with a single query specified as a string rather than a list
    patterns_to_match = []
    if type(patterns)==type("abc"):
        patterns_to_match.append(patterns)
        num_patterns=1
    else:
        patterns_to_match=patterns
    #create a pattern object for each query
    i=0
    for pattern in patterns_to_match:
        #if the pattern is too long simply abort, we don't want to deal
        #with overflow issues and there are better tools out there anyway
        #similarly if it less than 2 aa then dont proceed
        if not 1<len(pattern)<60:
            error = "Error: A query protein had a length outside the range supported by this program"
            raise Exception(error)
        matchers.append(compile(pattern,len(pattern),max_errors[i]))
        i+=1
    #open the specified database for reading    
    try:
        db_file = open(db_file_name, 'r')
    except IOError:
        error = "Error: Cannot open specified db file"
        raise Exception(error)     
    #storage for database information 
    rec = []
    header = ''
    #storage for information about each match
    #this allows for later calculation of p values and
    #output of this information to the user
    results_record_headers = []
    results_surrounding_text = []
    num_match_recs = 0
    #begin parsing database 
    for line in db_file:           
        #database body section
        if line[0] != '>':
            rec.append(line[:-1])               
        #database header
        else:  
            #attempt to match with previous body section
            record = ''.join(rec)
            begin_rec = 0
            for i in range(num_patterns):
                matches = agrepy(patterns_to_match[i], len(patterns_to_match[i]),\
                                 record,len(record),1,matchers[i])
                if matches:
                    if begin_rec==0:
                        begin_rec = 1
                        num_match_recs +=1
                    results_record_headers.append(split(header)[1])
                    pattern_matches = []
                    for elm in matches:
                        #for some large strings on 64 bit platforms sagrepy can return invalid
                        #match data, so remove this here rather than risk problems later
                        if elm[0]<elm[1]:
                            pattern_matches.append(record[elm[0]:elm[1]]) 
                    results_surrounding_text.append(pattern_matches)                  
            rec = []
            header = line
    #attempt to match the final body section
    record = ''.join(rec)
    begin_rec = 0
    for i in range(num_patterns):
        matches = agrepy(patterns_to_match[i],\
                    len(patterns_to_match[i]),record,len(record),1,matchers[i])
        if matches:
            if begin_rec==0:
                begin_rec = 1
                num_match_recs +=1
            results_record_headers.append(split(header)[1])
            pattern_matches = []
            for elm in matches:
                #for some large strings on 64 bit platforms sagrepy can return invalid
                #match data, so remove this here rather than risk problems later
                if elm[0]<elm[1]:
                    pattern_matches.append(record[elm[0]:elm[1]])  
            results_surrounding_text.append(pattern_matches)  
    #end parsing database
    #now calculate p values and output results to the user
    outputResults(patterns_to_match,num_match_recs,results_record_headers,results_surrounding_text)
    db_file.close()
    #Finally return result headers to the caller in case
    #they wish to make use of them in some other way
    return results_record_headers  
  
#------------------------------------------------------------------------------------------------------------------
#Statistical data on the background 'noise' of a typical FASTA database, computed from multiple runs of a
#a sample database (without filtering). These model the mean and standard deviation of the expected distribution #for each length query we could input for a match, from 1-79 inclusive. Note that the standard deviations are 
#still very wonky and are not currently used. TODO - fix generateBackgroundProfile, hopefully to generate these 
#on the fly for different databases. See generateBackgroundProfile for how these are computed.
means = [0.0576, 0.2060, 0.4227, 0.6167, 0.7813, 0.9330, 1.0194, 1.1055, 1.1532, 1.2232, 1.2708, 1.3107, 1.3779, 1.4190, 1.4751, 1.5233, 1.5680, 1.6106, 1.6532, 1.6971, 1.7290, 1.7653, 1.7979, 1.8408, 1.8838, 1.9175, 1.9429, 1.9744, 2.0082, 2.0253, 2.0481, 2.0724, 2.1001, 2.1122, 2.1365, 2.1729, 2.1979, 2.2104, 2.2427, 2.2739, 2.2839, 2.3123, 2.3111, 2.3388, 2.3891, 2.3998, 2.4093, 2.4425, 2.4651, 2.5089, 2.5061, 2.5610, 2.5678, 2.6084, 2.5905, 2.6712, 2.6958, 2.7556, 2.7995, 2.7929, 2.7638, 2.7955, 2.8336, 2.8287, 2.9153, 2.9848, 2.9831, 3.0588, 3.1133, 2.9859, 3.1840, 3.1121, 3.1644, 3.2074, 3.2510, 3.3560, 3.3798, 3.4434, 3.3666]

#---------------------------------------------------------------
# For all matches calculate p-values and print results to screen
#--------------------------------------------------------------- 
def outputResults(patterns_to_match,num_match_recs,results_record_headers,results_surrounding_text):
    total_pat = len(results_record_headers)
    if(num_match_recs>0):
        #score matches and generate p values
        rslt_scores = []
        pvalues = {}
        pat_nums = []
        max_elms = []
        cum_pat_num = 0
        i=0
        j=0
        while i<total_pat: 
            scr = 0
            pval = 1.0   
            pat_num = 0
            this_header=results_record_headers[i]
            while i<total_pat and results_record_headers[i]==this_header:
                #if there are multiple matches in this record for this pattern
                #then find the most likely and use its score
                max_score=0
                max_elm=0
                max_len=1
                elmnt = 0
                for elm in results_surrounding_text[i]:
                    len_sur_text = len(elm)
                    this_score = score(len_sur_text,patterns_to_match[pat_num],elm)   
                    if this_score>max_score:
                        max_score=this_score
                        max_len = len_sur_text
                        max_elm = elmnt
                    elmnt += 1
                max_elms.append(max_elm)
                #Assume values for rslt_score can be modelled by random variable X, 
                #having a normal distribution with mean M and standard deviation S. 
                #We get values for M and S from precomputed statistical data for 
                #this type of database and then we use them to normalize X 
                #to the standard normal distribution, by using the conversion 
                #Z = X - M / S, where Z is the random variable modelling the standard 
                #normal distribution. Thus, we can convert any result score value in
                #X to an equivalent value in Z, and this normalised value can be fed
                #into the Gaussian function, which uses the erfc function to find the
                #probability of X > rslt_score and X < -rslt_score, which is the
                #p value for this score. 
                #TODO use tables of means and standard deviations as follows:
                #gaussian((rslt_score-means[min_len-1])/stdevs[min_len-1])
                standardised_mean = max_score-means[max_len-1]
                if standardised_mean<0:
                    standardised_mean = 0
                pval *= gaussian(standardised_mean)
                scr += max_score
                pat_num+=1
                i+=1  
            cum_pat_num+=pat_num
            pat_nums.append(cum_pat_num) 
            rslt_scores.append(scr)            
            pvalues[j] =  pval
            j+=1
        #find ordering by pvals of results
        ordered = {}
        ordered = sorted(pvalues.items(),key=itemgetter(1))
        ordering = []
        cutoff = PVAL_THRESHOLD
        for elm in ordered:
            if elm[1]<=cutoff:
                ordering.append(elm[0])
            else:
                num_match_recs-=1
        #now print all this data, sorting according to p-values
        if(num_match_recs==1):
            print "One match was",
        else:
            print "In order of significance, these are the %i matches"%(num_match_recs)
        if len(patterns_to_match)>1:
            print "found containing the patterns in %s\n"%(patterns_to_match)
        else:
            print "found for the pattern %s\n"%(patterns_to_match)
        print "|  Record  |      In Text       |  Raw Score  |    P-Value    |" 
        print "|----------|--------------------|-------------|---------------|" 
        for i in range(num_match_recs):
            num = ordering[i]
            surd_text = []
            if num!=0:
                lower = pat_nums[num-1]
            else:
                lower=0
            header = results_record_headers[lower]
            #some databases enclose headers in braces so remove these
            header = header.strip('(').strip(')')
            for j in range(lower,pat_nums[num]):
                surd_text.append(results_surrounding_text[j][max_elms[j]])
                surd_text.append(" and ")
            sur_text = ''.join(surd_text[:-1])
            len_sur_text = len(sur_text)
            num_lines = (len_sur_text/20)+1
            middle_line = int(floor((num_lines+1)/2.0))-1
            for j in range(num_lines):
                if j==middle_line: 
                    print "|  %s "%(header),
                else:
                    print "|         ",
                if j==num_lines-1:
                    chars = len_sur_text%20
                    print "|%s"%(sur_text[j*20:j*20+chars]),
                    for k in range(19-chars):
                        print "",
                    print "|",
                else:        
                    print "|%s|"%(sur_text[j*20:(j+1)*20]),
                if j==middle_line:
                    print "     %.2i     |"%(rslt_scores[num]),
                    if pvalues[num]<9.999999e-99:
                        print " %e|"%(pvalues[num])
                    else:
                        print " %e |"%(pvalues[num])
                else:
                    print "            |               |"
            print "|----------|--------------------|-------------|---------------|"
    else:
        print "Sorry, could not find any results matching %s"%(patterns_to_match)
    print ""