print "->This is a quick demonstration of the Adrasteia program"
print "->Checking for string LILLFN in our sample database:"
print ""
import sys
sys.path.append('../Bin') 
from Adrasteia import *
direct_match("../Src/Testing/SampleDB.fasta","LILLFN",2)
print "->Now checking for database records that contain both LILLFN and PESKR"
print ""
direct_match("../Src/Testing/SampleDB.fasta",["LILLFN","PESKR"],2)
print "->Now a longer query, using preprocessing to speed up the search"
print "->This query has a match with one error and spanning two lines in"
print "->the database, and so it cannot be easily located by grep"
print ""
m=FastaDatabase("../Src/Testing/SampleDB.fasta")
m.match(["PAALTAVEMAGVKYLQVQHGSNVNIHR"],1)
print "->now testing for 3 fragments  with more errors using our preprocessed database"
print "->note that we specify a maximum of 3 errors, and one of the segments has 4"
print "->errors, thus discarding it from the result set"
print ""
m.match(["LDIKLIDYTM","IEERACEVNFISDKDLYVAAL","DEAMKPRSPSEYEDTSSPG"],3)
print "->Finally a more complicated batch query specifying grouping on some"
print "->of the elements and using different error rates for each"
print ""
m.match([["QMVEEADH","QTEETQKTVPEQ","ETQNTVEPEPTQE"],["TGTAHWLHNDGNT","AVIKR"],"QVNMIRHTIRPKGL"],[1,1,1,0,0,4])