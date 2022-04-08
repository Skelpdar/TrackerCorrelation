"""
Hacky script to correlate trigger scintillator, White Rabbit TDC and fiber tracker events

Author: Erik
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

#Various parameters
p = {
    "WRStartOfSpillPort": 1, #Physical port on the FMC TDC where the start of spill signal comes
    "WRTelescopePort": 2, #Physical port on the FMC TDC for telescope trigger signal
    #"TSFPGAFrequency":1.248623529*0.9995150078*0.9999702524*1.000012886*0.9999961275*0.9999993372/0.35528982*0.35528970*0.429724/0.42971903*0.428688011/0.428687812*0.444850191/0.444849690*1e-6/156., #Frequency of TS FPGA
    "TSFPGAFrequency":8.000046470962321e-09, #Frequency of TS FPGA
    "TimestampTolerance": 5e-7, #Difference (in seconds) for where two timestamps will be considered equal
    "FTToWRDistance": 0.0021-2.872e-6 #Time difference between telescope trigger signals in WR and FT
}

class TSSpill:
    def __init__(self, coarse):
        self.coarse = coarse #Rough time of start of spill spill in international atomic time
        self.ticks = [] #Clock ticks since start of spill for each event

    """
    Translates TS clock ticks to international atomic time modulo start of spill, using White Rabbit's measurements of the telescope scintillator trigger

    TODO: Make a robust algorithm to align the two datasets
    """
    def toIAT(self, WRSpills):

        spillCandidates = []
        for i in range(0,len(WRSpills)):
            if abs(WRSpills[i].startOfSpill-self.coarse) < 50:
                spillCandidates.append(i)

        for i in spillCandidates:
            WRSpill = WRSpills[i]

            for skip in [0,1,2,3,4]:
                    #bins = np.linspace(0,0.5,5000000)

                    #Keep track of how we have translated the WR timestamps
                    deltat = 0

                    tsnew = np.array(self.ticks)

                    tsnew = np.array(tsnew - min(tsnew))
                    deltat += min(WRSpill.timestamps)
                    timestampsnew = np.array(WRSpill.timestamps) - min(WRSpill.timestamps)

                    #Try skipping initial few WR timestamps
                    timestampsnew = timestampsnew[skip:]
                    
                    tsnew = (np.array(tsnew) - min(tsnew))
                    deltat += min(timestampsnew)
                    timestampsnew = np.array(timestampsnew) - min(timestampsnew)

                    #tsnew = tsnew*1.248623529*0.9995150078*0.9999702524*1.000012886*0.9999961275*0.9999993372/0.35528982*0.35528970
                    tsnew = tsnew*p["TSFPGAFrequency"]

                    #if(True):
                    #    bins = np.linspace(0.4,0.45,500000)

                    #    plt.hist(tsnew,histtype="step",bins=bins,density=True)
                    #    plt.hist(timestampsnew,bins=bins,histtype="step", density=True)
                    #    plt.show()

                    #Check whether we have a match
                    #closesttimestamp = []
                    matches = 0
                    misses = 0
                    for k in tsnew:
                        diff = np.sort(np.abs(timestampsnew-k))
                        if diff[1] < p["TimestampTolerance"]:
                            misses += 1
                            #closesttimestamp.append(-1)
                        elif diff[0] < p["TimestampTolerance"]:
                            matches += 1
                            d = np.abs(timestampsnew-k)
                            #closesttimestamp.append(float(timestampsnew[np.where(d == d.min())])+deltat+startofspill[i])
                        else:
                            misses += 1
                            #closesttimestamp.append(-1)

                    print("Efficiency: " + str(matches/len(tsnew)))

                    if matches/len(tsnew) < 0.5:
                        continue
                    matched = True
                    #results.append(matches/len(tsnew))
                    #print("Matched")

                    #print(tsnew)
                    #print(deltat)

                    return tsnew + deltat, i

                    break
        else:
            raise Exception("Could not correlate TS spill with WR")

"""
Takes a plaintext TS datafile and extracts all clock tick counts in event headers

I do not trust to never miss events. 
"""
def readTS(filename):
    file = open(filename)
    lines = []
    for l in file:
        lines.append(l)
    
    latestcoarse = 0

    spills = [TSSpill(0)]

    for i in range(0,len(lines)):
        #Find UTC time of data from packet header
        if "sw: got read " in lines[i]:
            value = int(lines[i][:-1].split()[-2])
            latestcoarse = value
        #Iterate through each event
        if "FFFFFFFFFFFFFFFF" in lines[i]:
            if len(lines[i+1]) == 17:
                value = str(bin(int(lines[i+1][8:],16))[2::].zfill(32))[::-1]
                #Extract to clock ticks in the event header
                #v = int(value[0:8][::-1] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2)
                v = int(value[0:8][::-1][6:] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2)

                if len(spills[-1].ticks) == 0:
                    spills[-1].ticks.append(v)
                    spills[-1].coarse = latestcoarse
                #Start a new spill
                elif v < spills[-1].ticks[-1]:
                    spills.append(TSSpill(latestcoarse))
                    spills[-1].ticks.append(v)
                else:
                    spills[-1].ticks.append(v)
            #The event is sometimes split up into two network packets, this skips the network packet header lines  
            elif len(lines[i+6]) == 17:
                value = str(bin(int(lines[i+6][8:],16))[2::].zfill(32))[::-1]
                #v = int(value[0:8][::-1] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2)
                v = int(value[0:8][::-1][6:] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2)
                if len(spills[-1].ticks) == 0:
                    spills[-1].ticks.append(v)
                    spills[-1].coarse = latestcoarse
                elif v < spills[-1].ticks[-1]:
                    spills.append(TSSpill(latestcoarse))
                    spills[-1].ticks.append(v)
                else:
                    spills[-1].ticks.append(v)

    return spills

class WRSpill:
    def __init__(self, startOfSpill):
        self.startOfSpill = startOfSpill
        self.timestamps = []

"""
Read a plaintext White Rabbit data file
"""
def readWR(filename):
    wrfile = open(filename)

    lines = []

    for l in wrfile:
        lines.append(l)

    spills = []

    for i in range(int(len(lines)/7)):
        if int(lines[i*7+2].split()[1]) == p["WRStartOfSpillPort"]:
            #Sometimes multiple start of spill signals arrive close together
            #so I don't start a new spill if they are less than 5 seconds apart
            if len(spills) == 0:
                spills.append(WRSpill(int(lines[i*7+4].split()[1])))
            elif int(lines[i*7+4].split()[1]) - spills[-1].startOfSpill > 5:
                spills.append(WRSpill(int(lines[i*7+4].split()[1])))
            
        if int(lines[i*7+2].split()[1]) == p["WRTelescopePort"]:
            #Subtract start of spill time, so we never have to think about floating point precision
            spills[-1].timestamps.append(float(lines[i*7+4].split()[1])+float(lines[i*7+5].split()[1])*8e-9+float(lines[i*7+6].split()[1])*8e-9/4096.-spills[-1].startOfSpill)

    return spills

"""
Data from a single fiber tracker (there are four of them)
"""
class FT:
    def __init__(self, offset):
        self.offset = offset
        self.timestamps = []
        self.data = []


    """
    Take a list of IAT timestamps
    """
    def separateTSEvents(self, tstimestamps, tsspill, wrStartOfSpill):
        ftnew = self.timestamps.copy()

        ftnew = np.array(ftnew) - wrStartOfSpill
        ftnew += self.offset

        #ftvnew += 0.0021-2.872e-6
        #ftvnew = ftvnew - min(ftvnew)

        ftnew = ftnew[ftnew > -10]
        ftnew = ftnew[ftnew < 10]

        #plt.vlines(ftvnew,0,1, label="FT Vertical")

        h = []

        matches = 0
        #misses = 0
        for i,k in enumerate(tstimestamps):
            hit = False

            #This is very slow
            diff = list(np.abs(ftnew-k))
            sort = np.sort(np.abs(ftnew-k))
            if sort[1] < p["TimestampTolerance"]:
                pass
            elif sort[0] < p["TimestampTolerance"]:
                hit = True
                matches += 1
                index = diff.index(sort[0])
                print(str(tsspill.ticks[i])+":"+("".join([str(s)+"," for s in self.data[index]])[:-1]))
                pass
            if hit == False:
                #print(tsspill.ticks[i])
                print(str(tsspill.ticks[i])+":"+"")

        #print("FTV efficiency: " + str(matches/len(TSTimestamps)))
        

def readFT(filename):
    file = open(filename)

    lines = []

    for l in file:
        lines.append(l)

    ft = FT(p["FTToWRDistance"])

    for i in range(len(lines)):
        s = lines[i].split()
        if len(s) == 3:
            ft.timestamps.append(float(s[1])+float(s[0])*8e-9)
            
            data = []
            for k in range(0,192):
                if s[2][k] == "1":
                    data.append(k)
            ft.data.append(data)

    return ft

if __name__ == '__main__':
    outputfilename = ""
    TSfilename = ""
    WRfilename = ""
    FTfilename = ""

    if len(sys.argv) == 2:
        print("Using default data paths")
        outputfilename = sys.argv[1]
        TSfilename = "data/positive4GevElectrons_Mar28_1455.txt"
        WRfilename = "data/March28_calibration_1.txt"
        FTfilename = "data/March28_1150_calibrationrun_1_51.txt"
    elif len(sus.argv) == 5:
        outputfilename = sys.argv[1]
        TSfilename = sys.argv[2]
        WRfilename = sys.argv[3]
        FTfilename = sys.argv[4]
    else:
        quit()

    print("Reading TS data")
    TSSpills = readTS(TSfilename)

    print("Reading WR data")
    WRSpills = readWR(WRfilename)

    print("Reading FT data")
    #Fiber tracker 51 is the downstream vertical tracker, 50 is the downstream horizontal
    FT_DownstreamVertical = readFT(FTfilename)

    #print(TSSpills[1].ticks)
     
    print("Correlating TS and WR and FT data")
    for i in range(0,len(TSSpills)):
        try:
            TSTimestamps, spillnumber = TSSpills[i].toIAT(WRSpills)
            print("Successfully correlated " + str(i))

            print("Correlating hits in fiber tracker with TS events")
            #This prints out a line for each TS event, with a comma separated list of hit positions in the fibers
            FT_DownstreamVertical.separateTSEvents(TSTimestamps, TSSpills[i], WRSpills[spillnumber].startOfSpill)
        except:
            print("Could not correlate " + str(i))
            print("Zeros for all TS events")
            for k in TSSpills[i].ticks:
                print(str(k)+":")
    
