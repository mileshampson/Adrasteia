Adrasteia 0.3.3
14.10.2006
Miles Hampson
------------------------------------------------------------------------
See the Legal subdirectory for copyright applying to this program.

For a quick demonstration run the demonstration.sh script, and to test
the program works run the testSuite.sh script.

On a Linux machines they will compile and run the program through a quick 
demonstration and test on your system. For the compilation to work you 
will need python 2.4 and swig 1.3 installed and callable on your machine.  
The compile script assumes that the python header files are in 
/usr/include/python2.4, if this is not the case on your machine please 
replace the 2 lines in the /src/testing/Compile.sh file of this package  
with the appropriate names.

------------------------------------------------------------------------
Supported Platforms and Architectures:

Has been tested to run on the Intel x86 and AMD64 architectures only
(although should run on the full range of 64 bit and PPC architectures).

The scripts have only been tested on linux, however the program has no
external dependencies apart from the standard python 2.4 and Swig 1.3
distributions,and so should compile on any platform with them installed. 

In order to compile for other platforms you will need to run the steps
int the /src/testing/Compile.sh for your OS, namely running swig on the
.i files to generate wrappers, then compiling your .c files and loading
the shared output files.

------------------------------------------------------------------------
Description of program:

Adrasteia is a simple program designed to provide an easy to use method
for searching for short protein strings in a FASTA database. It takes 
one or more peptide fragments, each specified as a string, and finds
areas in the database that closely resemble the query. How close can 
be specified by a parameter indicating the maximum difference between
the two. The results are ranked by how likely they are to be significant,
and displayed in order of significance. If more than one query is
specified and a database record has more than one inside it, its ranking
will be boosted. This means that you can input peptide fragments that
are possibly related, and have any records they are both present in 
come up top of the list of results.
-------------------------------------------------------------------------
How to use:

To quickly get going call direct_match(db_file_name,patterns,max_errors)
from the Adrasteia module. Specify the name of the database relative to
the current directory in db_file_name, and a single or list of patterns
to match in the database within the given max_error distance.

For larger databases it is recommenended that you first index the
database by creating an instance of FastaDatabase(db_file_name).
Once this is done you can call the much faster match() function.
-------------------------------------------------------------------------
match(patterns,max_errors)

Find matches in database within specified distance of the given pattern.
Parameters:
            
patterns    - a string (i.e. "MKFL")
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
            
max-errors  - an integer
            - a list (i.e [2,3,1])
            
            A single integer value will be applied to ALL patterns.
            
            If a seperate error rate is desired for each pattern then
            specify a list with 1 integer per pattern
------------------------------------------------------------------------
Indexing:

The first time a FastaDatabase instance is created for a database it
indexes the database and stores the results in .idx and .slice files
in the database directory. The process does about a Mb per minute for
every 256Mb of RAM in the computer, but will only need to be run once 
per database, after that the program loads the index files very quickly. 
The purpose of indexing is to store the distinguishing features of each 
record in a way that allows quick comparison to a query pattern, so that 
rather than having to run through the entire database with each pattern 
we wish to match we can jump directly to likely candidates. This works
especially well for multiple queries that are grouped together (see
above). For several queries matching will usually be an order of
magnitude quicker than a direct trawl through the database. The
trade off is the minute per Mb to first create the indices, and the
storage space for them (up to 2x the size of the database for long
word lengths). If this is not acceptable then either use the direct
match function which avoids indexing altogether, or adjust the
parameters of the indexing scheme (see DevelopersManual.pdf)

--------------------------------------------------------------------
Restrictions:

-The number of errors for a pattern must be between 1 and 8
-Use of match() requires:
    -patterns between (usual values) 2 and 60 aa in length
    -a length of at least twice the number of errors
--------------------------------------------------------------------