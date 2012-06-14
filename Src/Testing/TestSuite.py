print "->This test suite will test the Adrasteia program"
print "->against a very small sample database to check" 
print "->a set of input queries give expected results"
print "->This by no means is anywhere near complete coverage,"
print "->but should give confidence that data from most"
print "->possible input classes will be handled correctly"
print "->Due to the complexity of the results this testing"
print "->is not automated, you will have to eyeball the results"
print ""
import sys
sys.path.append('../Bin') 
from Adrasteia import *
m=FastaDatabase("../Src/Testing/SampleDB.fasta")

print "->1. Testing fragment directly matching middle of first record"
m.match(["DVNFISDKDLYVAALTNAD"],1)

print "->2. Testing fragment one deletion from middle of last record"
m.match(["EQNAKRIKDRTK"],1)

print "->3. Testing fragment one substitution from first characters of first record"
m.match(["MMFLILLFNILCLFP"],1)

print "->4. Testing fragment from first characters of last record"
print "with an extra first character added"
m.match(["DMENERAKQVYLAKLNE"],1)

print "->5. Testing fragment from last characters of first record,"
print "with last character deleted,"
m.match(["DEKERERFSI"],1)

print "->6. Testing fragment directly matching last characters of last record"
print "(which spans 2 lines)"
m.match(["LWTSDLEEGGK"],1)

print "->7. Testing small fragment that is not in the database"
m.match(["VQYFIKAGQ"],1)

print "->8. Testing large fragment that is not in the database"
m.match(["DEDTVDRLKNYNYRKTVRAKVCNRDRDYYYDSVIDKHEYYDDLVESTV"],1)

print "->9. Testing fragment that runs over 2 lines, with deletion"
print "on last character of line"
m.match(["QAFCNLKFRNVNRP"],1)

print "->10. Testing real fragment that has longest allowed length"
m.match(["GFLGGVAGVTSLLLMKASGTSMEEVRYWQYKWRLDRDENIQQAFKKLTEDENPELFKAHD"],1)

print "->11. Testing fragment that is 2 characters long, specified as string rather than list"
m.match("QA",1)

print "->12. Testing fragment that is one character long"
try:
    m.match(["Q"],1)
except Exception, inst:
    print inst
    print ""

print "->13. Testing empty fragment"
try:
    m.match([""],1)
except Exception, inst:
    print inst
    print ""

print "->14. Testing fragment with incorrect syntax"
try:
    m.match(1,["DVNFISDKDLYVAALTNAD"])
except Exception, inst:
    print inst
    print ""

print "->15. Testing fragment that matches within edit distance, but has invalid characters"
m.match(["MKF9ILLF"],1)

print "->16. Testing two fragments from the first record"
m.match(["QVLYESFNPLI","LYVAALTNADLNYTMVTPRP"],1)

print "->17. Testing two fragments from seperate records, one of which contains a deletion"
m.match(["SEAKTCENLVDTYRG","EPEQTEETQKTVEPEQTEE"],1)

print "->18. Testing two fragments that are not in database"
m.match(["VPDVEQTEVTQKTVVPETEE","LVVDQNEQ"],1)

print "->19. Testing two fragments only one of which is in database"
m.match(["DVNFISDKDLYVAALTNAD","LRGQVVPQF"],1)

print "->20. Testing two fragments that both run over 2 lines"
m.match(["SKSPRTASPTRRPSPKLP","KSPDEAMKRPRSPSEYED"],1)

print "->21. Testing two fragments that both run over 2 lines and share"
print "   the same middle line"
m.match(["VEMAGVKYLQVQHGSN","LYTGAIVTNNDGPYMAYVEVLGDPNLQFFIKSGDAWVTLSE"],1)

print "->22. Testing two fragments that are 1 character away from each other, and specified"
print "in the opposite order to which they appear in the database" 
m.match(["VVIKRASDRGFEWIAFKTNDNAITNLLAGRV","EVREGQVLMIPQN"],1)

print "->23. Testing two fragments that are adjacent, the first of which contains an insertion"
print "on the boundary between the two words"
m.match(["MENERAKQVYLAKLNEQAERYDEMVEAMKKVAALDVELTIM","EERNLLSVGYKNVIGARRASWRILSSIEQKEESKGN"],1)

print "->24. Testing two fragments that overlap"
m.match(["FYYK","KMKGDYFRYLAEFKSGA"],1)

print "->25. Testing three fragments on the same line"
m.match(["TTSESKVF","YYRYLAEFKIGDERK","TDLPPTHPIRLGLALNFS"],1)

print "->26. Testing a fragment in the database with an error rate of 0"
m.match(["AQLDVELTVEERNLVS"],0)

print "->27. Testing a fragment not in the database with an error rate of 0"
m.match(["AQLVELTEEEERLVS"],0)

print "->28. Testing a fragment with an error rate of 10"
m.match(["QAFDDAIAELDSLNEESYKDSTLIMQLLRDNLTLWT"],10)

print "->29. Testing fragment that is too long and with newlines in it"
try:
    m.match(["LLPWQKGQRSRPHHGHQQFQHQCDIQRLTASEPSRRVRS\nEAGVTEIWDHDTPEFRCAGFVAVRVVIQPGGLLLPSYSNAPYITFVEQGRGVQGVVVPGCPETFQSGSEFEYPRSQRDQRSRQSESGESSRGDQRSRQSESEESSRGDQRSRQSESEEFSRGDQHQKIFRIRDGDVIPSPAGVVQWTHNNGDNDLISITLYDANSFQNQLDEN\n\n\n\nVRNFFLAGQSKQSREDRRSQRQTREEGSDRQSRESQDDEALLEANILSGFEDEILQEIFRNVDQETISKLRGENDQRGFIVQARDLKLRVPEEYEEELQRERGDRKRGGSGRSNGLEQAFCNLKFRQNVNRPSRADVFNPRAGRINTVDSNNLPILEFIQLSAQHVVLYKNAILGPRWNLNAHSALYVTRGEGRVQVVGDEGRSVFDDNVQRGQILVVPQGFAVVLKAGREGLEWVELKNDDNAITSPIAGKTSVLRAIPVEVLANSYDISTKEAFRLKNGRQEVEVFRPFQSR"],99999)
except Exception, inst:
    print inst
    print ""

print "->30. Testing a fragment that is a header, with incorrectly specified error rates"
try:
    m.match([">14310_ARATH  P48347  14-3-3-like protein GF14"],[1,2])
except Exception, inst:
    print inst
    print ""

print "->31. Testing two fragments from seperate records that match with seperate error rates"
m.match(["YTMFHLADATYHECFKII","LQSIECQPQQSCTASL"],[1,2])

print "->32. Testing two fragments from seperate records one of whoose error rate doesnt match"
m.match(["YTMFHLADATYHECFKII","LQSIECQPQQSCTAS"],[2,1])

print "->33. Testing two grouped fragments from seperate records"
m.match([["YTMFHLADATYHECFKII","LQSIECQPQQSCTASL"]],[1,2])

print "->34. Testing two grouped fragments from the same record"
m.match([["IYDRNNGSIICLHLNYSPPSY","PEGPGASGLPPKAPGDK"]],[1,1])

print "->35. Testing four fragments grouped so they should match"
m.match([["IYDRNNGSIICLHLNYSPPSY","PEGPGASGLPPKAPGDK"],"YTMFHLADATYHECFKII","LQSIECQPQQSCTASL"],2)