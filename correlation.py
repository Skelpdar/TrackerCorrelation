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
    "HCalTimestampTolerance": 5e-6, #Difference (in seconds) for where two timestamps will be considered equal
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
        #for i in range(0,len(WRSpills)):
        #    if abs(WRSpills[i].startOfSpill-self.coarse) < 50:
        #        spillCandidates.append(i)

        spillCandidates = list(range(0,len(WRSpills)))
        print(spillCandidates)

        for i in spillCandidates:
            print(i)
            WRSpill = WRSpills[i]

            for skip in [0,1,2,3,4,5,6]:
                    try:
                        bins = np.linspace(0,0.5,5000000)

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

                        tsnew = tsnew*p["TSFPGAFrequency"]

                        #plt.hist(timestampsnew, bins=bins, histtype="step", linewidth=2)
                        #plt.hist(tsnew, bins=bins, histtype="step")
                        #plt.show()
                        
                        #Check whether we have a match
                        matches = 0
                        misses = 0
                        for k in tsnew:
                            diff = np.sort(np.abs(timestampsnew-k))
                            if diff[1] < p["TimestampTolerance"]:
                                misses += 1
                            elif diff[0] < p["TimestampTolerance"]:
                                matches += 1
                                d = np.abs(timestampsnew-k)
                            else:
                                misses += 1

                        print("Efficiency: " + str(matches/len(tsnew)))

                        if matches/len(tsnew) < 0.2:
                            continue
                        matched = True

                        return tsnew + deltat, i

                        break
                    except:
                        print("Exception?")
        else:
            raise Exception("Could not correlate TS spill with WR")

"""
Takes a plaintext TS datafile and extracts all clock tick counts in event headers

I do not trust this to never miss events. 
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
                v = int(value[0:8][::-1] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2)

                if v == 0:
                    value = str(bin(int(lines[i+2][8:],16))[2::].zfill(32))[::-1]
                    v = int(value[0:8][::-1] + value[8:16][::-1] + value[16:24][::-1]+ value[24:32][::-1],2) 

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
    def separateTSEvents(self, tstimestamps, tsspill, wrStartOfSpill, outputfile):
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
                outputfile.write(str(tsspill.ticks[i])+":"+("".join([str(s)+"," for s in self.data[index]])[:-1])+"\n")
                pass
            if hit == False:
                #print(tsspill.ticks[i])
                print(str(tsspill.ticks[i])+":"+"")
                outputfile.write(str(tsspill.ticks[i])+":"+"\n")

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

class HCalSpill:
    def __init__(self):
        self.timestamps = []

    def toIAT(self, WRSpills):
        spillCandidates = [0,1,2]
        #for i in range(0,len(WRSpills)):
        #    if abs(WRSpills[i].startOfSpill-self.coarse) < 50:
        #        spillCandidates.append(i)

        for i in spillCandidates:
            WRSpill = WRSpills[i]

            for skip in [0,1,2,3,4,5,6,7]:
                    bins = np.linspace(0,5,5000000)

                    #Keep track of how we have translated the WR timestamps
                    deltat = 0

                    tsnew = np.array(self.timestamps)

                    tsnew = np.array(tsnew - min(tsnew))
                    deltat += min(WRSpill.timestamps)
                    timestampsnew = np.array(WRSpill.timestamps) - min(WRSpill.timestamps)

                    #Try skipping initial few WR timestamps
                    timestampsnew = timestampsnew[skip:]
                    
                    tsnew = (np.array(tsnew) - min(tsnew))
                    deltat += min(timestampsnew)
                    timestampsnew = np.array(timestampsnew) - min(timestampsnew)

                    #tsnew = tsnew*p["TSFPGAFrequency"]
                   
                    plt.hist(tsnew, bins=bins, histtype="step",linewidth=2) 
                    plt.hist(timestampsnew, bins=bins, histtype="step") 
                    plt.show()
     
                    #Check whether we have a match
                    matches = 0
                    misses = 0
                    for k in tsnew:
                        diff = np.sort(np.abs(timestampsnew-k))
                        if diff[1] < p["HCalTimestampTolerance"]:
                            misses += 1
                        elif diff[0] < p["HCalTimestampTolerance"]:
                            matches += 1
                            d = np.abs(timestampsnew-k)
                        else:
                            misses += 1

                    print("Efficiency: " + str(matches/len(tsnew)))

                    if matches/len(tsnew) < 0.2:
                        continue
                    matched = True

                    return tsnew + deltat, i

                    break
        else:
            raise Exception("Could not correlate TS spill with WR")


def readHCal(filename):
    file = open(filename)

    spills = []

    oldtimestamp = 999999999
    for l in file:
        timestamp = float(l)/(5e6)
        if timestamp < oldtimestamp:
            spills.append(HCalSpill())
            spills[-1].timestamps.append(timestamp)
            oldtimestamp = timestamp
        else:
            spills[-1].timestamps.append(timestamp)
            oldtimestamp = timestamp
    
    return spills


if __name__ == '__main__':
    outputfilename = ""
    TSfilename = ""
    WRfilename = ""
    FTdvfilename = ""

    if len(sys.argv) == 2:
        print("Using default data paths")
        outputfilename = sys.argv[1]
        TSfilename = "data/ldmx_captan_out_17-04-2022_21-14-11__172.txt"
        WRfilename = "data/WR_out_172.txt"
        FTdvfilename = "data/DipClient_out_run_172_51.txt"
    elif len(sys.argv) == 5:
        outputfilename = sys.argv[1]
        TSfilename = sys.argv[2]
        WRfilename = sys.argv[3]
        FTdvfilename = sys.argv[4]
    else:
        print("Incorrect arguments")
        quit()

    print("Reading TS data")
    TSSpills = readTS(TSfilename)

    print(len(TSSpills))

    print("Reading WR data")
    WRSpills = readWR(WRfilename)
    print(len(WRSpills))

    #print("Reading HCal data")
    #HCalSpills = readHCal("data/hcal_172_timestamps.txt")

    #for i in range(0,len(HCalSpills)):
    #    iat = HCalSpills[i].toIAT([WRSpills[i]])

    print("Reading FT data")
    #Fiber tracker 51 is the downstream vertical tracker, 50 is the downstream horizontal
    #42 is upstream vertical, and 41 is upstream horizontal
    FT_DownstreamVertical = readFT(FTdvfilename)
    
    outputfile = open(outputfilename, "w")
   
    print(len(TSSpills[1].ticks))

    WRSpills = WRSpills[16:]    

     
    print("Correlating TS and WR and FT data")
    for i in range(0,len(TSSpills)):
        try:
            TSTimestamps, spillnumber = TSSpills[i].toIAT(WRSpills)
            print("Successfully correlated " + str(i))

            print("Correlating hits in fiber tracker with TS events")
            #This prints out a line for each TS event, with a comma separated list of hit positions in the fibers
            FT_DownstreamVertical.separateTSEvents(TSTimestamps, TSSpills[i], WRSpills[spillnumber].startOfSpill, outputfile)
        except:
            print("Could not correlate " + str(i))
            print("Zeros for all TS events")
            for k in TSSpills[i].ticks:
                outputfile.write(str(k)+":"+"\n")
    
