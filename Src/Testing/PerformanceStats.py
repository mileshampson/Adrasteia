print 'starting testing...'
import sys
sys.path.append('../Bin') 
from Adrasteia import *
sys.path.append('../Src/Statistics')
import Measurement
#Measurement.generateBackgroundProfile("DatabaseName",5,8)
#reference test, without preprocessing
Measurement.timeMultipleDirectMatch(1,5,24,26,"/home/tails/FullDB.fasta",1)
# test with preprocessing
Measurement.timeMultiplePreprocessAndMatch(1,1,24,26,"/home/tails/FullDB.fasta",1)
Measurement.timeMultipleDirectMatch(1,1,24,26,"/home/tails/FullDB.fasta",1)
#now run again, this time profiling the results of a filtered search
#import profile
#profile.run("db = FastaDatabase(\"DatabaseName\")")
#profile.run("db.match(\"PIRLGLALNFSVFYYEILNSPERACHLAKRAFDEAIAELDSLNEDSYKDSTLIMQLL\",3)")
#profile.run("direct_match(\"DatabaseName\",\"PIRLGLALNFSVFYYEILNSPER\",3)")
#import hotshot
#import hotshot.stats
#prof = hotshot.Profile("hotshot.txt")
#prof.runcall(FastaDatabase("DatabaseName"))
#prof.runcall(db.match("PIRLGLALNFSVFYYEILNSPERACHLAKRAFDEAIAELDSLNEDSYKDSTLIMQLL",3))
#prof.close()
#stats = hotshot.stats.load("hotshot.txt")
#stats.strip_dirs()
#stats.sort_stats('time', 'calls')
#stats.print_stats(20)