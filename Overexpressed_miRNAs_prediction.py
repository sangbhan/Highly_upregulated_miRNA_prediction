# This program finds overexpressed motifs and predicts overexpressed miRNAs for a given transcriptome data.

import sys, time

from scipy import stats

PSEUDO_COUNT   = 0.0000001
sNucList       = ["A", "C", "G", "T", "N"]
sHsChrList     = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]  
nSexChrDic     = {"X": 23, "Y": 24}
sStartCodon    = "ATG"
sStopCodons    = ["TAA", "TAG", "TGA"]
sComplementDic = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
fFoldChgCutOff = -0.50

# stores RefSeq information
class cRefSeq:

    # __init__ provides initial values for the instance variables of an object
    def __init__(self):
        
        self.sGeneSym      = ""
        self.sNMID         = ""
        self.nNMID         = 0
        self.nChrID        = 0
        self.sStrand       = ""
        self.nCDSStart     = 0
        self.nCDSEnd       = 0
        self.nExonCount    = 0
        self.nExonStarts   = []
        self.nExonEnds     = []
        self.nExonsSize    = 0
        self.n5UTRSize     = 0
        self.nORFSize      = 0
        self.n3UTRSize     = 0
        self.sExonsSeq     = ""
        self.s5UTRSeq      = ""
        self.sORFSeq       = ""
        self.s3UTRSeq      = ""
        self.n3UTRMotifDic = {}
        self.nORFMotifDic  = {}
        self.fFoldChg      = 0.0

    # end: __init__

    # parses each line in the RefFlat file 
    def ParseInfo(self, sInfoList):

        self.sGeneSym = sInfoList[0].upper()
        self.sNMID    = sInfoList[1]

        # removes "NM_" in the NMID and converts to int
        if "NM_" in self.sNMID:

            try:
                
                self.nNMID = int(self.sNMID.replace("NM_",""))

            except ValueError:

                print("\nThe RefFlat file contains an entry which has a non-integer value after NM_ in sNMID.\n")
            
                sys.exit()

        # end: if "NM_"

        # removes "chr" in the chromosome name and converts it to int
        if "chr" in sInfoList[2]:
            
            sChrID = sInfoList[2].replace("chr","")

            # converts sChrID to int
            if sChrID in sHsChrList:

                if sChrID in nSexChrDic.keys():

                    self.nChrID = nSexChrDic[sChrID]

                else:

                    self.nChrID = int(sChrID)

            # end: if sChrID
            
        # end: if "chr"

        self.sStrand = sInfoList[3]

        if (self.sStrand != "+") and (self.sStrand != "-"):

            print("\nThe RefFlat file contains an entry with wrong strand information.\n")

        try:

            self.nCDSStart  = int(sInfoList[6])
            self.nCDSEnd    = int(sInfoList[7])
            self.nExonCount = int(sInfoList[8])

        except ValueError:

            print("\nThe RefFlat file contains an entry with a non-integer value in nCdsStart, nCdsEnd, or nExonCount.\n")
            
            sys.exit()
        
        sExonStarts = sInfoList[9].split(",")[:-1]
        sExonEnds   = sInfoList[10].split(",")[:-1]

        if self.nExonCount != len(sExonStarts) or self.nExonCount != len(sExonEnds):

            print("\nThe RefFlat file contains an entry with wrong exon information.\n")

            sys.exit()

        # converts all strings in sExonStart and sExonEnds to int
        for nExonNum in range(self.nExonCount):

            try:

                self.nExonStarts.append(int(sExonStarts[nExonNum]))
                self.nExonEnds.append(int(sExonEnds[nExonNum]))

            except ValueError:

                print("\nThe RefFlat file contains a wrong entry with a non-integer number in nExonStarts or nExonEnds.\n")

                sys.exit()

        # end: for nExonNum

    # end: ParseInfo


    # stores exonic genome sequence, 5'UTR sequence, and 3'UTR sequence
    def StoreSeq(self, sSeq):

        # creates exonic genome sequence
        for nExonNum in range(self.nExonCount):

            self.sExonsSeq += sSeq[self.nExonStarts[nExonNum]:self.nExonEnds[nExonNum]]

            if self.nCDSStart in range(self.nExonStarts[nExonNum], self.nExonEnds[nExonNum]):

                nCDSStaExon = nExonNum

            if self.nCDSEnd - 1 in range(self.nExonStarts[nExonNum], self.nExonEnds[nExonNum]):

                nCDSEndExon = nExonNum

        # end: for nExonNum

        self.n5UTRSize = self.nCDSStart - self.nExonStarts[nCDSStaExon]
        self.n3UTRSize = self.nExonEnds[nCDSEndExon] - self.nCDSEnd

        for nExonNum in range(nCDSStaExon):

            self.n5UTRSize += self.nExonEnds[nExonNum] - self.nExonStarts[nExonNum]

        for nExonNum in range(nCDSEndExon + 1, self.nExonCount):

            self.n3UTRSize += self.nExonEnds[nExonNum] - self.nExonStarts[nExonNum]

        if self.sStrand == "-":

            self.sExonsSeq = Complement(self.sExonsSeq)
            self.n5UTRSize, self.n3UTRSize = self.n3UTRSize, self.n5UTRSize

        self.nExonsSize = len(self.sExonsSeq)
        self.nORFSize   = self.nExonsSize - self.n5UTRSize - self.n3UTRSize

        if self.nORFSize < 0:

            print("\nThe RefFlat file contains an entry with wrong ORF information.\n")
        
        self.s5UTRSeq = self.sExonsSeq[:self.n5UTRSize]
        self.s3UTRSeq = self.sExonsSeq[-self.n3UTRSize:]

        if self.n3UTRSize == 0:

            self.sORFSeq = self.sExonsSeq[self.n5UTRSize:]

        else:

            self.sORFSeq = self.sExonsSeq[self.n5UTRSize:-self.n3UTRSize]

    # end: StoreSeq


    # counts 3'UTR 7mer motifs
    def Count3UTRMotif(self, cMotifDic):

        # counts 3'UTR 7mer motifs
        for nPos in range(self.n3UTRSize - 6):

            try:

                self.n3UTRMotifDic[self.s3UTRSeq[nPos:nPos + 7]] += 1

            except KeyError:

                self.n3UTRMotifDic[self.s3UTRSeq[nPos:nPos + 7]]  = 1

                try:

                    cMotifDic[self.s3UTRSeq[nPos:nPos + 7]].append(self)

                except KeyError:

                    cMotifDic[self.s3UTRSeq[nPos:nPos + 7]]  = [self]

        # end: for nPos

        return cMotifDic

    # end: Count3UTRMotif


    # counts ORF 7mer motifs
    def CountORFMotif(self, cMotifDic):

        # counts ORF 7mer motifs
        for nPos in range(self.nORFSize - 6):

            try:

                self.nORFMotifDic[self.sORFSeq[nPos:nPos + 7]] += 1

            except KeyError:

                self.nORFMotifDic[self.sORFSeq[nPos:nPos + 7]]  = 1

                try:

                    cMotifDic[self.sORFSeq[nPos:nPos + 7]].append(self)

                except KeyError:

                    cMotifDic[self.sORFSeq[nPos:nPos + 7]]  = [self]

        # end: for nPos

        return cMotifDic

    # end: CountORFMotif


    # adds 1 to N1 or N3 and counts the number of highly downregulated genes and the number of not highly downregulated genes 
    def AddN1N3(self, c3UTRMotifDic, cORFMotifDic, nDown, nNotDown):

        # adds 1 to N1 and counts the number of highly downregulated genes
        if self.fFoldChg < fFoldChgCutOff:

            nDown += 1
            
            for sMotif in self.n3UTRMotifDic.keys():

                c3UTRMotifDic[sMotif].AddN1()

            for sMotif in self.nORFMotifDic.keys():

                cORFMotifDic[sMotif].AddN1()

        # end: if self.fFoldChg

        # adds 1 to N3 and counts the number of not highly downregulated genes
        else:

            nNotDown += 1
            
            for sMotif in self.n3UTRMotifDic.keys():

                c3UTRMotifDic[sMotif].AddN3()

            for sMotif in self.nORFMotifDic.keys():

                cORFMotifDic[sMotif].AddN3()

        # end: else

        return nDown, nNotDown

    # end: AddN1N3


    # sets n3UTRMotifDic = {}, nORFMotifDic = {}, and fFoldChg = 0.0
    def ResetMotifInfo(self):

        self.n3UTRMotifDic = {}
        self.nORFMotifDic  = {}
        self.fFoldChg      = 0.0

    # end: ResetMotifInfo

                    
    # sets fFoldChg
    def SetFoldChg(self, fFoldChg):

        self.fFoldChg = fFoldChg

    # end: SetFoldChg


    # gets sGeneSym
    def GetGeneSym(self):

        return self.sGeneSym

    # end: GetGeneSym


    # gets sNMID
    def GetsNMID(self):

        return self.sNMID

    # end: GetsNMID


    # gets nNMID
    def GetnNMID(self):

        return self.nNMID

    # end: GetnNMID


    # gets nChrID
    def GetChrID(self):

        return self.nChrID

    # end: GetChrID


    # gets sORFSeq
    def GetORFSeq(self):

        return self.sORFSeq

    # end: GetORFSeq


    # gets n5UTRSize
    def Get5UTRSize(self):

        return self.n5UTRSize

    # end: Get5UTRSize


    # gets nORFSize
    def GetORFSize(self):

        return self.nORFSize

    # end: GetORFSize


    # gets n3UTRSize
    def Get3UTRSize(self):

        return self.n3UTRSize

    # end: Get3UTRSize


    # gets fFoldChg
    def GetFoldChg(self):

        return self.fFoldChg

    # end: GetFoldChg

# end: cRefSeq


# stores each motif's information
class cORFMotif:

    # __init__ provides initial values for the instance variables of an object
    def __init__(self, sNewMotif):

        self.sMotif      = sNewMotif
        self.fPValue     = 0.0
        self.fBonferroni = 0.0
        self.nN1         = 0
        self.nN2         = 0
        self.nN3         = 0
        self.nN4         = 0
        self.fRelRisk    = 0.0

    # end: __init__


    # calculates nN2 and nN4
    def CalN2N4(self, nDown, nNotDown):

        self.nN2 = nDown - self.nN1
        self.nN4 = nNotDown - self.nN3

    # end: CalN2N4


    # calculates fRelRisk
    def CalRelRisk(self):

        if (self.nN1 + self.nN2 != 0) and (self.nN3 + self.nN4 != 0) and (self.nN3 != 0):

            try:

                fA = (self.nN1 + PSEUDO_COUNT) / (self.nN1 + self.nN2 + PSEUDO_COUNT)
                fB = (self.nN3 + PSEUDO_COUNT) / (self.nN3 + self.nN4 + PSEUDO_COUNT)

                self.fRelRisk = fA / fB

                if self.fRelRisk < 0.0:

                    print("\nRelative risk is less than 0.0.\n")

                    sys.exit()

            except ZeroDivisionError:

                print("\nZeroDivisionError occurs while calculating relative risk.\n")

                sys.exit()

    # end: CalRelRisk


    # calculates P Value            
    def CalPValue(self):

        fOddsRatio, self.fPValue = stats.fisher_exact([[self.nN1, self.nN2], [self.nN3, self.nN4]])

        if self.fPValue < 0.0 or self.fPValue > 1.0:

            print("\nP Value is less than 0.0 or more than 1.0.\n")
            sys.exit()

    # end: CalPValue


    # corrects P value by using the bonferonni correction
    def Bonferroni(self, nTestNum):

        self.fBonferroni = self.fPValue * nTestNum

    # end: Bonferroni


    # Predicts the overexpressed miRNAs
    def FindMiR(self, cM8Dic, cA1Dic):

        s7merM8MiR = ""
        s7merA1MiR = ""

        try:
            
            # creates s7merM8MiR
            for cMiR in cM8Dic[self.sMotif]:

                if s7merM8MiR == "":

                    s7merM8MiR = cMiR.GetMiRSym()

                else:

                    s7merM8MiR += ", " + cMiR.GetMiRSym()

            # end: for cMiR

            try:

                # creates s7merA1MiR
                for cMiR in cA1Dic[self.sMotif]:

                    if s7merA1MiR == "":

                        s7merA1MiR = cMiR.GetMiRSym()

                    else:

                        s7merA1MiR += ", " + cMiR.GetMiRSym()

                # end: for cMiR

                return "7mer-m8: " + s7merM8MiR + ", 7mer-A1: " + s7merA1MiR                

            except KeyError:

                return "7mer-m8: " + s7merM8MiR

        except KeyError:

            try:

                # creates s7merA1MiR
                for cMiR in cA1Dic[self.sMotif]:

                    if s7merA1MiR == "":

                        s7merA1MiR = cMiR.GetMiRSym()

                    else:

                        s7merA1MiR += ", " + cMiR.GetMiRSym()

                # end: for cMiR

                return "7mer-A1: " + s7merA1MiR

            except KeyError:

                return ""
            
    # end: FindMiR


    # adds 1 to nN1
    def AddN1(self):

        self.nN1 += 1

    # end: AddN1


    # adds 1 to nN2
    def AddN2(self):

        self.nN2 += 1

    # end: AddN2


    # adds 1 to nN3
    def AddN3(self):

        self.nN3 += 1

    # end: AddN3


    # adds 1 to nN4
    def AddN4(self):

        self.nN4 += 1

    # end: AddN4


    # gets sMotif
    def GetMotif(self):

        return self.sMotif

    # end: GetMotif


    # gets fBonferroni
    def GetBonferroni(self):

        return self.fBonferroni

    # end: GetBonferroni


    # gets nN1
    def GetN1(self):

        return self.nN1

    # end: GetN1


    # gets nN2
    def GetN2(self):

        return self.nN2

    # end: GetN2


    # gets nN3
    def GetN3(self):

        return self.nN3

    # end: GetN3


    # gets nN4
    def GetN4(self):

        return self.nN4

    # end: GetN4


    # gets fRelRisk
    def GetRelRisk(self):

        return self.fRelRisk

    # end: GetRelRisk


# end: cORFMotif


# stores each motif's information
class c3UTRMotif:

    # __init__ provides initial values for the instance variables of an object
    def __init__(self, sNewMotif):

        self.sMotif      = sNewMotif
        self.fPValue     = 0.0
        self.fBonferroni = 0.0
        self.nN1         = 0
        self.nN2         = 0
        self.nN3         = 0
        self.nN4         = 0
        self.fRelRisk    = 0.0

    # end: __init__


    # calculates nN2 and nN4
    def CalN2N4(self, nDown, nNotDown):

        self.nN2 = nDown - self.nN1
        self.nN4 = nNotDown - self.nN3

    # end: CalN2N4


    # calculates fRelRisk
    def CalRelRisk(self):

        if (self.nN1 + self.nN2 != 0) and (self.nN3 + self.nN4 != 0) and (self.nN3 != 0):

            try:

                fA = (self.nN1 + PSEUDO_COUNT) / (self.nN1 + self.nN2 + PSEUDO_COUNT)
                fB = (self.nN3 + PSEUDO_COUNT) / (self.nN3 + self.nN4 + PSEUDO_COUNT)

                self.fRelRisk = fA / fB

                if self.fRelRisk < 0.0:

                    print("\nRelative risk is less than 0.0.\n")

                    sys.exit()

            except ZeroDivisionError:

                print("\nZeroDivisionError occurs while calculating relative risk.\n")

                sys.exit()

    # end: CalRelRisk


    # calculates P Value            
    def CalPValue(self):

        fOddsRatio, self.fPValue = stats.fisher_exact([[self.nN1, self.nN2], [self.nN3, self.nN4]])

        if self.fPValue < 0.0 or self.fPValue > 1.0:

            print("\nP Value is less than 0.0 or more than 1.0.\n")
            sys.exit()

    # end: CalPValue


    # corrects P value by using the bonferonni correction
    def Bonferroni(self, nTestNum):

        self.fBonferroni = self.fPValue * nTestNum

    # end: Bonferroni


    # Predicts the overexpressed miRNAs
    def FindMiR(self, cM8Dic, cA1Dic):

        s7merM8MiR = ""
        s7merA1MiR = ""

        try:
            
            # creates s7merM8MiR
            for cMiR in cM8Dic[self.sMotif]:

                if s7merM8MiR == "":

                    s7merM8MiR = cMiR.GetMiRSym()

                else:

                    s7merM8MiR += ", " + cMiR.GetMiRSym()

            # end: for cMiR

            try:

                # creates s7merA1MiR
                for cMiR in cA1Dic[self.sMotif]:

                    if s7merA1MiR == "":

                        s7merA1MiR = cMiR.GetMiRSym()

                    else:

                        s7merA1MiR += ", " + cMiR.GetMiRSym()

                # end: for cMiR

                return "7mer-m8: " + s7merM8MiR + ", 7mer-A1: " + s7merA1MiR                

            except KeyError:

                return "7mer-m8: " + s7merM8MiR

        except KeyError:

            try:

                # creates s7merA1MiR
                for cMiR in cA1Dic[self.sMotif]:

                    if s7merA1MiR == "":

                        s7merA1MiR = cMiR.GetMiRSym()

                    else:

                        s7merA1MiR += ", " + cMiR.GetMiRSym()

                # end: for cMiR

                return "7mer-A1: " + s7merA1MiR

            except KeyError:

                return ""
            
    # end: FindMiR


    # adds 1 to nN1
    def AddN1(self):

        self.nN1 += 1

    # end: AddN1


    # adds 1 to nN2
    def AddN2(self):

        self.nN2 += 1

    # end: AddN2


    # adds 1 to nN3
    def AddN3(self):

        self.nN3 += 1

    # end: AddN3


    # adds 1 to nN4
    def AddN4(self):

        self.nN4 += 1

    # end: AddN4


    # gets sMotif
    def GetMotif(self):

        return self.sMotif

    # end: GetMotif


    # gets fBonferroni
    def GetBonferroni(self):

        return self.fBonferroni

    # end: GetBonferroni


    # gets nN1
    def GetN1(self):

        return self.nN1

    # end: GetN1


    # gets nN2
    def GetN2(self):

        return self.nN2

    # end: GetN2


    # gets nN3
    def GetN3(self):

        return self.nN3

    # end: GetN3


    # gets nN4
    def GetN4(self):

        return self.nN4

    # end: GetN4


    # gets fRelRisk
    def GetRelRisk(self):

        return self.fRelRisk

    # end: GetRelRisk


# end: c3UTRMotif


# stores each miRNA's information
class cMiRBase:

    # __init__ provides initial values for the instance variables of an object
    def __init__(self):

        self.sMiRSym     = ""
        self.sMiRSeq     = ""
        self.sComplement = ""
        self.s7merM8     = ""
        self.s7merA1     = ""

    # end: __init__


    # parses each line in the mature miRNA fasta file 
    def ParseInfo(self, sMiRSym, sMiRSeq):

        self.sMiRSym     = sMiRSym
        self.sMiRSeq     = sMiRSeq
        self.s7merM8 = Complement(sMiRSeq[1:8].replace("U","T"))
        self.s7merA1 = Complement(sMiRSeq[1:7].replace("U","T")) + "A"

    # end: ParseInfo


    # gets sMiRSym
    def GetMiRSym(self):

        return self.sMiRSym

    # end: GetMiRSym


    # gets s7merM8
    def Get7merM8(self):

        return self.s7merM8

    # end: Get7merM8


    def Get7merA1(self):

        return self.s7merA1

    # end: Get7merA1

# end: cMiRBase  


# removes non-NM sequences
def RemoveNonNM(cRefList):

    cNMList = []

    # keeps NM sequences
    for cRef in cRefList:

        if cRef.GetnNMID() != 0:

            cNMList.append(cRef)

    # end: for cRef

    return cNMList

# end: RemoveNonNM


# keeps NM sequences aligned only on chr 1-22, X, or Y
def KeepAutoXY(cRefList): 

    cAutoXYList = []
    
    # keeps NM sequences aligned only on chr 1-22, X, or Y
    for cRef in cRefList:

        if cRef.GetChrID() != 0:

            cAutoXYList.append(cRef)

    # end: for cRef

    return cAutoXYList

# end: KeepAutoXY


# removes NM sequences that have multiple entries in the RefFlat file
def RemoveRedun(cRefList):

    nNMIDDic    = {}
    cUniqueList = []

    # counts the numbers of sNMID in cRefList
    for cRef in cRefList:

        try:

            nNMIDDic[cRef.GetsNMID()] += 1

        except KeyError:

            nNMIDDic[cRef.GetsNMID()]  = 1

    # end: for cRef
        
    # keeps NM sequences that have a unique entry in the RefFlat file
    for cRef in cRefList:

        if nNMIDDic[cRef.GetsNMID()] == 1:

            cUniqueList.append(cRef)

    # end: for cRef

    return cUniqueList

# end: RemoveRedun


# checks the validity of sSeq
def ChrSeqCheck(sSeq, nChrNum):

    # checks the validity of sSeq
    for nPos in range(len(sSeq)):

        if sSeq[nPos] not in sNucList:

            print("\nThere is an abnormal character in the human chromosome", sHsChrList[nChrNum - 1], "at nucleotide", nPos + 1, "\n")

            sys.exit()

    # end: for nPos

# end: ChrSeqCheck


# converts a DNA sequence into its complement
def Complement(sSeq):

    sComplement = ""

    for nNucPos in range(len(sSeq)):

        sComplement += sComplementDic[sSeq[-nNucPos-1]]

    return sComplement

# end: Complement


# removes NM sequences that have wrong ORFs
def RemoveWrong(cRefList):

    cORFOKList = []

    # appends NM sequences that have wrong ORFS: no start codon, no stop codon, ORF size of 3N+1 or 3N+2, or internal stop codons to cORFOKList
    for cRef in cRefList:

        # appends NM sequences to cORFOKList if they have wrong ORFs
        if (cRef.GetORFSeq()[:3] == sStartCodon) and (cRef.GetORFSeq()[-3:] in sStopCodons) and (cRef.GetORFSize()%3 == 0):

            cORFOKList.append(cRef)

            # removes NM sequence that have internal stop codons
            for nCodonPos in range(1, int(cRef.GetORFSize()/3)-1):

                if cRef.GetORFSeq()[3*nCodonPos:3*nCodonPos+3] in sStopCodons:

                    del cORFOKList[-1]

                    break

            # end: for nCodonPos

        # end: if cRef.GetORFSeq()[:3], cRef.GetORFSeq()[-3:], and cRef.GetORFSize()%3

    # end: for cRef

    return cORFOKList

# end: RemoveWrong


# For each gene that has multiple isoforms, selects a representative isoform
def SelectRepNM(cRefList):

    cGeneSymDic = {}
    cRepRefList = []

    # creates cGeneSymDic with keys from cRef.GetGeneSym() and values set to [cRefs]
    for cRef in cRefList:

        try:

            cGeneSymDic[cRef.GetGeneSym()] += [cRef]

        except KeyError:

            cGeneSymDic[cRef.GetGeneSym()]  = [cRef]

    # end: for cRef

    # appends a representative isoform to cRepRefList
    for sGeneSym in cGeneSymDic.keys():

        cRepRef = cGeneSymDic[sGeneSym][0]

        # selects a representative isoform
        for nRefNum in range(1, len(cGeneSymDic[sGeneSym])):

            if cGeneSymDic[sGeneSym][nRefNum].GetnNMID() < cRepRef.GetnNMID():

                cRepRef = cGeneSymDic[sGeneSym][nRefNum]

        # end: for nRefNum

        cRepRefList.append(cRepRef)

    # end: for sGeneSym

    return cRepRefList

# end: SelectRepNM


# assigns the provided log2(fold-change) values to each cRef 
def FoldChgRef(cRefList, nDataNum):

    try:

        InFile = open("Mission5_Dataset" + str(nDataNum + 1) + ".txt", "r")

    except IOError:

        print("\nMission5_Dataset" + str(nDataNum + 1) + ".txt: No such file or directory\n")

        sys.exit()

    fFoldChgDic = {}

    for sReadLine in InFile.readlines():

        sInfoList = sReadLine.split()
        
        fFoldChgDic[sInfoList[0].upper()] = float(sInfoList[1])

    InFile.close()

    cFoldChgList = []

    # assigns the provided log2(fold-change) values to each cRef 
    for cRef in cRefList:

        if cRef.GetGeneSym() in fFoldChgDic.keys():

            cFoldChgList.append(cRef)

            cRef.SetFoldChg(fFoldChgDic[cRef.GetGeneSym()])

    # end: for cRef
    
    return cFoldChgList

# end: FoldChgRef


# uses the provided log2(fold-change) values and discovers the motifs that are significantly enriched in the highly downregulated genes after miR-1 transfection into HeLa cell.

def FindDownMotifs(cRefList, nDataNum, c7merM8Dic, c7merA1Dic):

    cFoldChgRef = FoldChgRef(cRefList, nDataNum)

    c3UTRMotifDic = {}
    cORFMotifDic  = {}

    for cRef in cFoldChgRef:

        c3UTRMotifDic = cRef.Count3UTRMotif(c3UTRMotifDic)
        cORFMotifDic  = cRef.CountORFMotif(cORFMotifDic)

    c3UTRObjectDic = {}
    cORFObjectDic  = {}

    for sMotif in c3UTRMotifDic.keys():

        cNewMotif = c3UTRMotif(sMotif)
        c3UTRObjectDic[sMotif] = cNewMotif

    for sMotif in cORFMotifDic.keys():

        cNewMotif = cORFMotif(sMotif)
        cORFObjectDic[sMotif] = cNewMotif

    nDownRef    = 0 
    nNotDownRef = 0

    for cRef in cFoldChgRef:

        nDownRef, nNotDownRef = cRef.AddN1N3(c3UTRObjectDic, cORFObjectDic, nDownRef, nNotDownRef)

        cRef.ResetMotifInfo()

    c3UTREnriched = []
    cORFEnriched  = []
    
    n3UTRNumOfTests = 0
    nORFNumOfTests  = 0

    # calculates N2, N4, RelRisk, and PValue of cMotifObject in c3UTRObjectDic.values()
    for cMotifObject in c3UTRObjectDic.values():

        cMotifObject.CalN2N4(nDownRef, nNotDownRef)
        cMotifObject.CalRelRisk()

        # picks and counts the motifs that have A/B > 1.0
        if cMotifObject.GetRelRisk() > 1.0:

            n3UTRNumOfTests += 1
            
            cMotifObject.CalPValue()
            c3UTREnriched.append(cMotifObject)

        # end: cMotifObject.GetRelRisk()

    for cMotifObject in c3UTREnriched:

        cMotifObject.Bonferroni(n3UTRNumOfTests)

    c3UTREnriched.sort(key = c3UTRMotif.GetBonferroni)

    # end: for cMotifObject

    # calculates N2, N4, RelRisk, and PValue of cMotifObject in cORFObjectDic.values()
    for cMotifObject in cORFObjectDic.values():

        cMotifObject.CalN2N4(nDownRef, nNotDownRef)
        cMotifObject.CalRelRisk()

        # picks and counts the motifs that have A/B > 1.0
        if cMotifObject.GetRelRisk() > 1.0:

            nORFNumOfTests += 1
            
            cMotifObject.CalPValue()
            cORFEnriched.append(cMotifObject)

        # end: cMotifObject.GetRelRisk()

    # end: for cMotifObject

    for cMotifObject in cORFEnriched:

        cMotifObject.Bonferroni(nORFNumOfTests)

    cORFEnriched.sort(key = cORFMotif.GetBonferroni)

    print("\nDataSet" + str(nDataNum + 1))
    print("3'UTR")

    # reports results for 3'UTR 7mers
    for nTopNumber in range(5):

        try:

            print(c3UTREnriched[nTopNumber].GetMotif(), c3UTREnriched[nTopNumber].GetBonferroni(), c3UTREnriched[nTopNumber].GetN1(), c3UTREnriched[nTopNumber].GetN2(), c3UTREnriched[nTopNumber].GetN3(), c3UTREnriched[nTopNumber].GetN4(), c3UTREnriched[nTopNumber].GetRelRisk(), c3UTREnriched[nTopNumber].FindMiR(c7merM8Dic, c7merA1Dic), sep = "\t", end = "\n")

        except IndexError:

            print("\nLess than 5 enriched motifs in the highly downregulated genes")
            break

    # end: for nTopNumber

    print("ORF")

    # reports results for ORF 7mers
    for nTopNumber in range(5):

        try:

            print(cORFEnriched[nTopNumber].GetMotif(), cORFEnriched[nTopNumber].GetBonferroni(), cORFEnriched[nTopNumber].GetN1(), cORFEnriched[nTopNumber].GetN2(), cORFEnriched[nTopNumber].GetN3(), cORFEnriched[nTopNumber].GetN4(), cORFEnriched[nTopNumber].GetRelRisk(), cORFEnriched[nTopNumber].FindMiR(c7merM8Dic, c7merA1Dic), sep = "\t", end = "\n")

        except IndexError:

            print("\nLess than 5 enriched motifs in the highly downregulated genes")
            break

    # end: for nTopNumber

# end: FindDownMotifs


# main function
def main():

    print("The program starts at", time.ctime())

    cRefSeqList = []

    try:

        InFile = open("refFlat.txt", "r")

    except IOError:

        print("\nrefFlat.txt: No such file or directory\n")

        sys.exit()
        
    # parses each line in the RefFlat file 
    for sReadLine in InFile.readlines():

        sInfoList = sReadLine.split()

        if len(sInfoList) != 11:

            print("\nThe RefFlat file contains an entry with more than or less than 11 values\n")
            
            sys.exit()

        cReadRefSeq = cRefSeq()
        cReadRefSeq.ParseInfo(sInfoList)
            
        cRefSeqList.append(cReadRefSeq)

    # end: for sReadLine

    InFile.close()

    cNMRefSeq  = RemoveNonNM(cRefSeqList)
    cAutoXYRef = KeepAutoXY(cNMRefSeq)
    cUniqueRef = RemoveRedun(cAutoXYRef)

    cChrRefDic = {}

    # creates cChrRefDic with keys from cRef.GetChrID() and values set to [cRefs]
    for cRef in cUniqueRef:

        try:

            cChrRefDic[cRef.GetChrID()] += [cRef]

        except KeyError:

            cChrRefDic[cRef.GetChrID()]  = [cRef]

    # end: for cRef

    # reads a chromosome sequence file and stores exonic genome sequence, 5'UTR sequence, and 3'UTR sequence of cRef
    for nChrNum in sorted(cChrRefDic.keys()):

        try:

            InFile = open("chr" + sHsChrList[nChrNum - 1] + ".fa", "r")

            InFile.readline()
            sChrSeq = InFile.read()

            InFile.close()

            sChrSeqPro = sChrSeq.replace("\n", "").upper()

            ChrSeqCheck(sChrSeqPro, nChrNum)

            for cRef in cChrRefDic[nChrNum]:

                cRef.StoreSeq(sChrSeqPro)

        except IOError:

            print("\nchr" + sHsChrList[nChrNum - 1] + ".fa: No such file or directory\n")

            sys.exit()

    # end: for nChrNum

    cORFOKRef   = RemoveWrong(cUniqueRef)
    cRepRefSeq  = SelectRepNM(cORFOKRef)

    cMiRList = []

    try:

        InFile = open("mature.fa", "r")

    except IOError:

        print("\nmature.fa: No such file or directory\n")

        sys.exit()

    sMiRInfoList = InFile.readlines()

    # parses each line in the mature miRNA fasta file
    for nMiRNum in range(int(len(sMiRInfoList) / 2)):

        # parses each line in the mature miRNA fasta file
        if "hsa" in sMiRInfoList[2 * nMiRNum].split()[0]:

            cReadMiRBase = cMiRBase()
            cReadMiRBase.ParseInfo(sMiRInfoList[2 * nMiRNum].split()[4], sMiRInfoList[2 * nMiRNum + 1].replace("\n", "").upper())

            cMiRList.append(cReadMiRBase)

        # end: if "hsa"
            
    # end: for nMiRNum

    c7merM8Dic = {}
    c7merA1Dic = {}

    # creates c7merM8Dic and c7merA1Dic
    for cMiR in cMiRList:

        try:
            
            c7merM8Dic[cMiR.Get7merM8()].append(cMiR)

        except KeyError:

            c7merM8Dic[cMiR.Get7merM8()] = [cMiR]

        try:
            
            c7merA1Dic[cMiR.Get7merA1()].append(cMiR)

        except KeyError:

            c7merA1Dic[cMiR.Get7merA1()] = [cMiR]

    # end: for cMiR

    for nDataSetNum in range(3):

        FindDownMotifs(cRepRefSeq, nDataSetNum, c7merM8Dic, c7merA1Dic)

    print("\nThe program ends at", time.ctime())
    
# end: main


main()
