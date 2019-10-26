import numpy as np
import sys
import math
import os
import signal
import time

# class to read input card
class inReader:
    def __init__(self,inFile):
        with open(inFile,'r') as f:
            content = f.readlines()
        f.close()
        self.data = content
        
        # data structure to be added
        self.numPin = 0
        self.name = []
        self.ifctype = []
        self.fArray = []

        self.sInput = ''
        self.fInput = []
        self.ifcPre = ''
        
        self.sPath = ''
        self.sOmp = ''
        self.sMpi = ''

        self.fPath = ''

        self.maxIter = 0
        self.mCouple = 0
        self.epsi = 0.0

        self.segCon = []
        self.segArray = []
        self.axial = []
        self.angular = []
        self.hollow = []

        self.energy = []

    ########################################################################## 
    # run the reading process
    ########################################################################## 
    def run(self):
        inReader.addName(self)
        
        inReader.addSinput(self)
        inReader.addFinput(self)
    
        inReader.addIFC(self)
        inReader.addIfcPre(self)
    
        inReader.addFpath(self)
    
        inReader.addSpath(self)
        inReader.addOmp(self)
        inReader.addMpi(self)

        inReader.add_mCouple(self)
        inReader.addMaxIter(self)
        inReader.addEpsi(self)
    
        inReader.addSeg(self)
        inReader.addSegArray(self)
    
        inReader.addAxial(self)
        inReader.addAngular(self)
        inReader.addHollow(self)
        inReader.addEnergy(self)
    
    ########################################################################## 
    # read input files
    ########################################################################## 
    def addSinput(self):
        for line in self.data:
            if 'sInput' in line and not line.startswith('#'):
                sInput = line.split()[-1]
        self.sInput = self.sInput + sInput

    def addFinput(self):
        for line in self.data:
            if 'fInput' in line and not line.startswith('#'):
                fInput = line.split()[1:]
        self.fInput = self.fInput + fInput

    def addIfcPre(self):
        for line in self.data:
            if 'ifcPre' in line and not line.startswith('#'):
                ifcPre = line.split()[-1]
        self.ifcPre = self.ifcPre + ifcPre

    ########################################################################## 
    # read serpent infomation
    ########################################################################## 
    # serpent run path
    def addSpath(self):
        for line in self.data:
            if 'sPath' in line and not line.startswith('#'):
                sPath = line.split()[-1]
        self.sPath = self.sPath + sPath
    # number of omp
    def addOmp(self):
        for line in self.data:
            if 'omp' in line and not line.startswith('#'):
                nOmp = int(line.split()[-1])
        if  nOmp > 1:
            sOmp = '-omp '+str(nOmp)
        else:
            sOmp = ''

        self.sOmp = self.sOmp + sOmp 

    # number of mpi
    def addMpi(self):
        for line in self.data:
            if 'mpi' in line and not line.startswith('#'):
                nMpi = int(line.split()[-1])
        if  nMpi > 1:
            sMpi = '-mpi '+str(nMpi)
        else:
            sMpi = ''

        self.sMpi = self.sMpi + sMpi

    ########################################################################## 
    # read FEMAXI data
    ########################################################################## 
    def addFpath(self):
        for line in self.data:
            if 'fPath' in line and not line.startswith('#'):
                fPath = line.split()[-1]
        self.fPath = self.fPath + fPath

    ########################################################################## 
    # read coupling interface control data
    ##########################################################################
    # read coupled method
    def add_mCouple(self):
        for line in self.data:
            if 'mCouple' in line and not line.startswith('#'):
                mCouple = int(line.split()[-1])
        self.mCouple = self.mCouple + mCouple
    def addMaxIter(self):
        for line in self.data:
            if 'maxIter' in line and not line.startswith('#'):
                maxIter = int(line.split()[-1])
        self.maxIter = self.maxIter + maxIter

    def addEpsi(self):
        for line in self.data:
            if 'epsi' in line and not line.startswith('#'):
                epsi = float(line.split()[-1])
        self.epsi = self.epsi + epsi


    ########################################################################## 
    # add ifc file name
    ########################################################################## 
    def addFileIFC(self,fArray):
        self.fArray = self.fArray + fArray

    ########################################################################## 
    # read pin geometry
    ########################################################################## 
    # read pin name and number
    def addName(self):
        name = []
        for line in self.data:
            if 'cPin' in line and not line.startswith('#'):
                for i in range(len(line.split())):
                    if i != 0:
                        name.append(line.split()[i])
        self.name = self.name + name
        self.numPin = self.numPin + len(name) 
    # read pin ifc type
    def addIFC(self):
        ifc = []
        for line in self.data:
            if 'ifctype' in line and not line.startswith('#'):
                for i in range(self.numPin+1):
                     if i != 0:
                         ifc.append(line.split()[i])
        self.ifctype = self.ifctype + ifc
    
    # read segmentation information
    def addSeg(self):
        segCon = []
        for line in self.data:
            if 'segPin' in line and not line.startswith('#'):
                for i in range(self.numPin+1):
                     if i != 0:
                         segCon.append(line.split()[i])
        self.segCon = self.segCon + segCon

    # read name of segments for each pin
    def addSegArray(self):
        segArray = []
        for i in range(len(self.data)):
            line = self.data[i]
            if 'segArray' in line and not line.startswith('#'):
                for j in range(1,self.numPin+1):
                    data = self.data[i+j].split()
                    segArray.append(data)
        self.segArray = self.segArray + segArray

    # read axial data for each pin
    def addAxial(self):
        axial = []
        for i in range(len(self.data)):
            line = self.data[i]
            if 'axial' in line and not line.startswith('#'):
                for j in range(1,self.numPin+1):
                    data = self.data[i+j].split()
                    axial.append(data)
        self.axial = self.axial + axial

    # read angular data for each pin
    def addAngular(self):
        angular = []
        for i in range(len(self.data)):
            line = self.data[i]
            if 'angular' in line and not line.startswith('#'):
                for j in range(1,self.numPin+1):
                    data = self.data[i+j].split()
                    angular.append(data)
        self.angular = self.angular + angular

    # read radial data for each pin
    def addHollow(self):
        hollow = []
        for i in range(len(self.data)):
            line = self.data[i]
            if 'hollow' in line and not line.startswith('#'):
                for j in range(1,self.numPin+1):
                    hollow.append(float(self.data[i+j]))
        self.hollow = self.hollow + hollow

    # read fast neutron energy data for each pin
    def addEnergy(self):
        energy = []
        for i in range(len(self.data)):
            line = self.data[i]
            if 'energy' in line and not line.startswith('#'):
                for j in range(1,self.numPin+1):
                    data = self.data[i+j].split()
                    energy.append(data)
        self.energy = self.energy + energy


              
########################################################################## 
# class for pin data
########################################################################## 
class pinData:
    # initialize pin data
    def __init__(self,name, ifctype, segPin):
        self.name = name
        self.ifctype = ifctype
        self.segPin = segPin 
       
        # Geometry group, values to be added later  
        self.nAxial = 0
        self.axial_range = []
        self.axial = []

        self.nCrack = 0
        self.angular_range = []
        self.angular = []

        self.nRadial = 0
        self.hollow = 0.0
        self.radial_range = []
        self.angular_range = []
        self.radial = []
       
        self.energy_range = []
        # Axial segmentation group, values to be added later  
        self.segName = []

    ########################################################################## 
    # pin geometry structure
    ########################################################################## 
    # add number of axial nodes per pin
    def addNaxial(self,nAxial):
        self.nAxial = self.nAxial + int(nAxial)
    # add axial range  data
    def addAxial_Range(self,axial_range):
        self.axial_range = self.axial_range + axial_range

    # add number of crackers per pin
    def addNcrack(self,nCrack):
        self.nCrack = self.nCrack + int(nCrack)
    # add agular range data
    def addAngular_Range(self,angular_range):
        self.angular_range = self.angular_range + angular_range

    # add central hole radius data
    def addHollow(self,hollow):
        self.hollow = self.hollow + hollow
    # add number of radial nodes per pin
    def addNradial(self,nRadial):
        self.nRadial = self.nRadial + int(nRadial)
    # add radial range  data
    def addRadial_Range(self,radial_range):
        self.radial_range = self.radial_range + radial_range

    # add fast neutron energy
    def addEnergy_Range(self,energy_range):
        self.energy_range = self.energy_range + energy_range


    def addAxial(self,axial):
        self.axial = self.axial + axial
    def addAngualr(self,angular):
        self.angular = self.angular + angular
    def addRadial(self,radial):
        self.radial = self.radial + radial

    ########################################################################## 
    # pin geometry structure
    ########################################################################## 
    # judge wether axial segmentation is actived
    def segPin_bool(self):
        if self.segPin == 1:
            return True
        else:
            return False
    def addSeg(self,segName):
        self.segName = self.segName + segName

    ########################################################################## 
    # FEMAXI fuel temperature array
    ##########################################################################
    def addFuelT(self,fuelT):
        self.fuelT = fuelT

    ########################################################################## 
    # read pin data from the class
    ########################################################################## 
    # get axial nodes numbers 
    def get_n_axial(self):
        return self.nAxial
    # get angular ring crack numbers 
    def get_n_crack(self):
        return self.nCrack
    # get angular ring crack numbers 
    def get_n_crack(self):
        return self.nCrack

    # get axial range
    def get_axial_range(self):
        return self.axial_range
    # get angular range
    def get_angular_range(self):
        return self.angular_range
    # get radial range
    def get_radial_range(self):
        return self.radial_range

    # get radial centre hole diameter
    def get_hollow(self):
        return self.hollow
    # get segmenation pin name
    def get_segName(self):
        return self.segName



########################################################################## 
# class to read FEMAXI output data
########################################################################## 
class FEMAXI_reader:
    def __init__(self):
    #    self.fuelT = {}
        self.n_axial = 0
        self.n_radial = 0
        self.time = []
        self.suffix = []
        self.rmax = 0.0
        self.rmax_gas = 0.0
        self.rmax_clad = 0.0

    # function for find time value in FEMAXI output
    def find_between_(s, first, last):
        try:
            start = s.index( first ) + len( first )
            end = s.index( last, start )
            return s[start:end]
        except ValueError:
            return ""

    # function for reading data from FEMAXI output
    def read(self, infile):
        # change infile to the right format
        infile = infile+'.out'
        # time array
        time = []
        suffix = []
        # pop out timepoint index array
        popidx = []
        # read steady state operation time from input
        with open(infile,'r') as fin:
            for line in fin:
                if '1 ' and 'HISTORY DATA (1)' in line:
                   # print (line)
                    for i in range(8):
                        line = next(fin)
                    # read steady state time
                    while True:
                        line = fin.readline()
                        if not line.strip():
#                            print ('pass')
                            pass
                        elif line.split()[-1] == '0':
#                            print ('iteration step',line.split()[0])
                            newline = line.replace(':',' ')
                            timeburn_str = newline.split()
                        elif line.split()[-1] == '2': 
                            break
                    timeburn = float(timeburn_str[1])*3600+float(timeburn_str[2])*60+float(timeburn_str[3])+float(timeburn_str[4])/1000
                    break    
        print ('operational time (s): ', timeburn)
#        print ('operational time: ', timeburn)
      #  print ('operational time: ', timeburn_str)
      #  sys.exit()
    
                
        # read fuel tempearture array from input
        for line in open(infile,'r'):
        	# read time
        	if 'TIME (H:M:S:MS)' in line:
                    time_str = FEMAXI_reader.find_between_(line,')','|')         
                    time_str = time_str.split(':')
                # print (time_str)
        
        	    # convert time to s
                    time_curr = float(time_str[0])*3600+float(time_str[1])*60+float(time_str[2])+float(time_str[3])/1000 - timeburn
        	    #	print (time_curr)
                    suffix_curr = str('{0:.2f}'.format(time_curr)).replace('.','d') # convert total time to transient time
        	#	print (suffix_curr)
                    time.append('%.2f'%float(time_curr))
                    suffix.append(suffix_curr)
    #    print (time)
    #    print (suffix)
        for timepoint in time:
            if float(timepoint) < 0:
                popidx.append(time.index(timepoint))      
    #    print ('pop out list is:', popidx)  
        # read pin geometry data 
        with open(infile,'r') as fin:
            for line in fin:
                if 'PELLET SPECIFICATIONS ' in line:
                    line = next(fin)
                    line = next(fin)
                    line = next(fin)
                    dPellet = float(line.split()[2])
                    rmax = dPellet/2
                elif 'CLAD. INSIDE DIAMETER (CM)' in line:
                    dGas = float(line.split()[-1])
                    dClad = float(next(fin).split()[-1])
                    rmax_gas = dGas/2
                    rmax_clad = dClad/2
        fin.close()
    
        # read feul temperature nodes
        with open(infile,'r') as fin:
            for line in fin:
                if 'TEMPERATURE DISTRIBUTION IN THE FUEL (DEG.C) ' in line:
                    # print (line)
                    line = next(fin)
                    line = next(fin)
                    line = next(fin)
                    # print (line)
                    data = line.split()
                    n_axial = int(data[0])
                    n_radial = len(data) - 1
                    break
            # add gas and cladding node to radcal array
            n_radial = n_radial + 2 # 1 node for gas and 1 node for cladding
        fin.close()
        # print(n_axial,n_radial)
        
        # create temperature array
        
        # number of time point
        n_old = len(time)
        
        # creat dictionary for fuel temperature
        fuelT = {}
        
        with open(infile,'r') as fin:
            n = 0
            # read fuel temperature into dictionary
            for line in fin:
                if 'TEMPERATURE DISTRIBUTION IN THE FUEL (DEG.C) ' in line:
                    data = []
                    T_idx = 'fuelT'+str(n)
                    line = next(fin)
                    line = next(fin)
                    for i in range(n_axial):
                        data.append((next(fin).split()))
                    data = np.array(data, dtype = float)
                    data = data + 273.15 # convert C to K
                    fuelT[T_idx] = data[:,1:n_radial+1]
                    # print (fuelT[T_idx])
                    n = n+1
    
            # read cladding and gas temperature into dictionary
            fin.seek(0)
            n = 0
            for line in fin:
                if 'PC      PS      CI      CO ' in line:
                    Tgas = []
                    Tclad = []
                    T_idx = 'fuelT'+str(n)
    
                    for i in range(n_axial):
                        data = next(fin)
                        data = data.replace('*',' ')
                        Tgas.append((float(data.split()[8])+float(data.split()[9]))/2)
                        Tclad.append((float(data.split()[9])+float(data.split()[10]))/2)
    
                    Tgas = np.array(Tgas, dtype = float)
                    Tclad = np.array(Tclad, dtype = float)
                    
                    Tgas = Tgas + 273.15 # convert C to K
                    Tclad = Tclad + 273.15 # convert C to K
    
                    # append gas and clad temperature to the fuel temperature array
                    fuelT[T_idx] = np.c_[fuelT[T_idx],Tgas]
                    fuelT[T_idx] = np.c_[fuelT[T_idx],Tclad]
                    n = n+1
                  #  fuelT[T_idx] = np.append(fuelT[T_idx],Tclad,axis=1)
       # sys.exit()
        fin.close()
    
        # remove element from burnup mode
        time = [i for j, i in enumerate(time) if j not in popidx]
        # convert time array to float
        time = [float(i) for i in time]
    #    print(type(time[2]))	
    #    print ('time is ',time)
        suffix = [i for j, i in enumerate(suffix) if j not in popidx]
    	
    #    for key in fuelT.keys():
    #        print ('old fuelT key are: ',key)
    #    print ('length of old fuelT key: ',len(fuelT))
        for i in popidx:
            del_idx = 'fuelT'+str(i)
            del fuelT[del_idx]
    #    print ('fuelT is ',fuelT)
    #    print ('time is ',time)
    #    print ('suffix is ', suffix)
        
        # replace the old key by renamed new key
        n_rm = len(popidx)
        Tidx_new = range(n_old - n_rm -1)
    #    print (len(fuelT))
    #    print (Tidx_new)
        for Tidx_new in range(n_old - n_rm):
            key_new = 'fuelT'+str(Tidx_new)
            key_old = 'fuelT'+str(Tidx_new+n_rm)
       #     print(key_new)
       #     print(key_old)
            fuelT[key_new] = fuelT.pop(key_old)
    #    for key in fuelT.keys():
    #        print ('new fuelT keys are: ',key)
    #    print ('length of new fuelT key: ',len(fuelT))
        self.fuelT = fuelT
        self.n_axial = self.n_axial + n_axial
        self.n_radial = self.n_radial + n_radial
        self.time = self.time + time
        self.suffix = self.suffix + suffix
        self.rmax = self.rmax + rmax
        self.rmax_gas = self.rmax_gas + rmax_gas
        self.rmax_clad = self.rmax_clad + rmax_clad

    # extact data from class
    def get_fuelT(self):
        return self.fuelT
    def get_n_axial(self):
        return self.n_axial
    def get_n_radial(self):
        return self.n_radial
    def get_time(self):
        return self.time
    def get_suffix(self):
        return self.suffix
    def get_rmax(self):
        return self.rmax
    def get_rmax_gas(self):
        return self.rmax_gas
    def get_rmax_clad(self):
        return self.rmax_clad

########################################################################## 
# class to generate IFC file for each pin 
########################################################################## 
class IFCgen():
    def __init__(self):
        self.time = []
        self.suffix = []
        pass

    # function for calculate the segment height
    def seg_Axial(zmin, n_axial, infile):
            infile = str(infile) + '.out'
            z = np.zeros(n_axial+1)
            z[0] = zmin
    #        z[n_axial] = zmax
            with open(str(infile), 'r') as inf:
                for line in inf:
                    if '*INPUT DATA' in line:
                        line = next(inf)
                        line = next(inf)
                        line = next(inf)
                        for i in range(n_axial):
                            dL = next(inf).split()[-1]
                            z[i+1] = z[i] + float(dL)
           # print (z)
            return z
    
    # function for calculate the angel degree of each pellet
    def pellet_Angle(n_Crack,angular):
            a = np.zeros(n_Crack+1)
            a[0] = angular[0]
            a[n_Crack] = angular[1]
            if n_Crack > 1:
                dA = (amax - amin)/n_Crack
                for i in range(n_Crack - 1):
                    a[i+1] = a[i] + dA 
            return a
        
     # function for calculate out radius of each fuel ring
    def ring_Radius(n_radial, rmin, rmax, rmax_gas,rmax_clad):
        r = np.zeros(n_radial+1)
       # print (rmin)
        r[0] = rmin
    #	r[n_radial] = rmax
        r[-3] = rmax        # outer fuel pellet boundary
        PI = 3.141592653	# constant PI
        S_pin = PI*(rmax**2 - rmin**2)
        S_par = S_pin/(n_radial-2) # fuel pellet rings
        for i in range(n_radial-1-2):
            r[i+1] =  '%.4f'%math.sqrt(S_par/PI + r[i]**2)
       # print (r)
        r[-2] = rmax_gas    # inner cladding boundary
        r[-1] = rmax_clad   # outer cladding boundary
       # print (r)
       # sys.exit()
        return r[1:n_radial+1+2]

    # calculate pin radial range
    def r_range(r,hollow):
        r_range = []
        r_range.append(hollow)
        r_range.append(r[-1])

        return r_range
    
    
    def dataGen(pinList,inData):
        fDataList = []
        # FAMXI data list
        for i in range(len(pinList)):
            Name = 'fData'+str(i+1)
            fDataList.append(Name)
   
        # looping to generate dict for pin overall, axial and radical data 
        for i in range(len(pinList)):
            # read FEMAXI data
            fDataList[i] = FEMAXI_reader()
            FEMAXI_reader.read(fDataList[i],inData.fInput[i])

            axial = IFCgen.seg_Axial(pinList[i].axial_range[0], fDataList[i].n_axial, inData.fInput[i])
        #    print (axial)
            angular = IFCgen.pellet_Angle(pinList[i].nCrack, pinList[i].angular_range)
        #    print (angular)
            radial = IFCgen.ring_Radius(fDataList[i].n_radial,pinList[i].hollow,fDataList[i].rmax,fDataList[i].rmax_gas,fDataList[i].rmax_clad)
        #    print (radial)

            r_range = IFCgen.r_range(radial,pinList[i].hollow)  # calculate radial range

            pinList[i].addNaxial(fDataList[i].n_axial)
            pinList[i].addNradial(fDataList[i].n_radial)
        #    print (pinList[i].nRadial)
        #    print (pinList[i].nAxial)
            pinList[i].addRadial_Range(r_range)
       #     print (pinList[i].radial_range)

            pinList[i].addAxial(list(axial))
            pinList[i].addAngualr(list(angular))
            pinList[i].addRadial(list(radial))
      #      print (pinList[i].axial)
      #      print (pinList[i].angular)
      #      print (pinList[i].radial)

            pinList[i].addFuelT(fDataList[i].fuelT)
   
#        self.time = self.time + fDataList[0].time
#        self.suffix = self.suffix + fDataList[0].suffix
        time = fDataList[0].time
        suffix = fDataList[0].suffix

        return time, suffix

       # for timestep indd range(len(self.suffix)):
    def writeData(suf,pinList,numPin,ifcPre):
        # write data into Serpent ifc file of each time step
        # print (suffix)
        # print (timestep)
        # define the outfile name
        outfile_step = ifcPre+'_' + suf
        # fuel array index
        fuelTkey = 'fuelT0'#+str(timestep)
        #    print ('current key is: ',fuelTkey)
        # print (outfile_step)
        with open(outfile_step, 'w+') as fout:
           # print('start writing: ',outfile_step)
           # fout.write( str(ifctype) +' ' + str(outfile) + ' ' + str(len(pinname))+ '\n' )
            fout.write(str(pinList[0].ifctype) +' ' + 'OUT_'+str(outfile_step) + ' ' + str(numPin)+ '\n' ) # try both format
            # write each pin data to file
            for i in range(numPin):
               # print (pinname.split())
               # print (n)
               # fout.write(pinname[n].split()[1] + '\n') # pin universe name
                if pinList[i].segPin:
                    fout.write ('-'+str(len(pinList[i].segName))+' ')
                    for j in range(len(pinList[i].segName)):
                        fout.write(str(pinList[i].segName[j])+' ')
                    fout.write('\n')
                else:
                    fout.write(pinname[n] + '\n') # pin universe name
                # segmating data of the pin
    #            fout.write(' '.join(map(str,pin[pinname[n]])) + '\n')
    #            fout.write(' '.join(map(str,pin[pinname[n]])) + ' ' + str(emin) +' '+ str(emax) + '\n') # whether this line is manditory is not clear
                fout.write(str(pinList[i].nAxial)+' ')
                fout.write(str(pinList[i].axial_range[0])+' ')
                fout.write(str(pinList[i].axial_range[1])+' ')
    
                fout.write(str(pinList[i].nCrack)+' ')
                fout.write(str(pinList[i].angular_range[0])+' ')
                fout.write(str(pinList[i].angular_range[1])+' ')
    
                fout.write(str(pinList[i].nRadial)+' ')
                fout.write(str(pinList[i].radial_range[0])+' ')
                fout.write(str(pinList[i].radial_range[1])+' ')
                fout.write('\n')

                fout.write(str(pinList[i].nAxial)+' ')
                fout.write(str(pinList[i].axial_range[0])+' ')
                fout.write(str(pinList[i].axial_range[1])+' ')
    
                fout.write(str(pinList[i].nCrack)+' ')
                fout.write(str(pinList[i].angular_range[0])+' ')
                fout.write(str(pinList[i].angular_range[1])+' ')
    
                fout.write(str(pinList[i].nRadial)+' ')
                fout.write(str(pinList[i].radial_range[0])+' ')
                fout.write(str(pinList[i].radial_range[1])+' ')

                fout.write(str(pinList[i].energy_range[0])+' ')
                fout.write(str(pinList[i].energy_range[1])+' ')
                fout.write('\n')

                fout.write(str(pinList[i].nAxial)+'\n')
        
                # looping segment from bottom to top, 0 angel to 360 angel, inner ring to out ring
                for j in range(pinList[i].nAxial):
                    fout.write(str('%.4f'%pinList[i].axial[j]) + ' ' + str('%.4f'%pinList[i].axial[j+1]) + ' ' + str(pinList[i].nCrack) + '\n') # axial segment dat
                    for k in range(pinList[i].nCrack):
                        fout.write(str(pinList[i].angular[k]) + ' ' + str(pinList[i].angular[k+1]) + ' ' + str(pinList[i].nRadial) + '\n') # angel segment data
                        for l in range(pinList[i].nRadial):
                            fout.write(str('%.4f'%pinList[i].radial[l]) + ' ' + str('%.4f'%pinList[i].radial[l]) + ' ' + str('%.2f'%pinList[i].fuelT[fuelTkey][pinList[i].nAxial-j-1][l]) + '\n') # radical segment data
    
        fout.close()

            
    def gen(pinList,inData):
        time, suffix = IFCgen.dataGen(pinList,inData)
        numPin = inData.numPin
        ifcPre = inData.ifcPre
        fArray = [] # array for ifc file names
        # current scheme suits for one-way coupling, to be improved 
        for timestep in range(len(suffix)):
            suf = suffix[timestep]
            IFCgen.writeData(suf,pinList,numPin,ifcPre)
            ifcFile = ifcPre+'_'+suf
            fArray.append(ifcFile)
        inData.addFileIFC(fArray)

# update serpent input file with coupled interface and signal control
class sAdd:
    def __init__(self,sInput,ifcPre,maxIter):
        ifcFile = ifcPre+'_'+'file'
        with open (sInput,'r') as f:
            inFile = f.readlines()
        f.close()
        sInput = sInput+'_run'
        with open (sInput,'w') as f:
            for line in inFile:
                f.write(line)
        f.close()

        self.sInput = sInput
        self.inFile = inFile
        self.ifcFile = ifcFile
        self.inSignal = 'com.in'
        self.outSignal = 'com.out'
        self.maxIter = maxIter

    def addSignalMode(self):
        with open (self.sInput,'a') as f:
            f.write('set comfile'+' ')
            f.write(self.inSignal+' ')
            f.write(self.outSignal+' ')
            f.write('\n')

    def addIFC(self):
        with open (self.sInput,'a') as f:
            f.write('ifc ')
            f.write(str(self.ifcFile))
            f.write('\n')
    def addMaxIter(self):
        with open(self.sInput,'a') as f:
            f.write('set ccmaxiter'+' ')
            f.write(str(self.maxIter))
            f.write('\n')
    def run(self):
        sAdd.addSignalMode(self)
        sAdd.addIFC(self)
        sAdd.addMaxIter(self)


# run serpent and FEMAXI
class excute:
    def sRun(sPath,sInput,omp,mpi):
        sInput = sInput+'_run'
        sRun = sPath+' '+sInput+' '+omp+' '+mpi
        run = 'nohup ' +sRun+' &> out.log & echo $! > pid.log'  # run background, log file stored in out.log, pid stored in pid.log
        os.system(run)
    def fRun(self):
        pass

# control of the coupled simulation
class control:
    def __init__(self,fArray,sInput):
        self.simulating = True
        self.iterating = True
        self.sleeping = True
        self.inSignal = 'com.in'
        self.outSignal = 'com.out'

        self.sInput = sInput
        self.fArray = fArray
        self.fTime = []
        self.sTime = []
        self.ifcFile = ''
        self.conv_P = ''
        self.conv_T = ''

        self.Power = []
        self.p_Error = []

    # steady state coupling
    def iterateSS(self,pinList,inData):
        control.convFile(self)
        epsi = inData.epsi
        while self.iterating:
            ###################
            # Wait for signal #
            ###################
            control.waitSignal(self)
            # Check if iteration has finished 
            if not self.iterating:
                break

            ###########################
            # Check convergency       #
            ###########################
            converge = control.tConverge(self,epsi)
 
            if not converge:
                ###########################
                # Read power distribution #
                ###########################
                control.pRead(self,self.fArray[0])
                ###########################
                # Run FEMAXI              #
                ###########################
                # generate new ifc data
                IFCgen.gen(pinList,inData)
                ###########################
                # Update Fuel Temperature #
                ###########################
                control.writeIFC(self.fArray[0])
                time.sleep(2)
                # Signal Serpent (SIGUSR1) #
                control.inSignal_1(self)
            else:
                # Signal Serpent (SIGUSR2) #
                control.inSignal_2(self)
                # Signal Serpent (SIGTERM) #
                control.outSignal_term(self)

    # one-way transient coupling
    def iterateOneway(self):
        control.convFile(self)
    #    epsi = inData.epsi
        # generate serpent time steps
        control.sTimeGen(self)
        control.fTimeGen(self)
        control.reDup(self)

      #  print (self.fTime)
      #  print (self.sTime)

        step = 0
        while self.simulating:
            self.iterating = True
            step += 1
            if not self.simulating:
                break
            while self.iterating:
               # print (self.iterating,step)
                ###################
                # Wait for signal #
                ###################
                control.waitSignal(self)
                # Check if iteration has finished 
                if not self.iterating:
                    break
                if self.sTime[step] in self.fTime:
                    print ('update ifc',self.sTime[step])
                    idx = self.fTime.index(self.sTime[step])
                #    print (idx)
                    ###########################
                    # Update Fuel Temperature #
                    ###########################
                    control.writeIFC(self.fArray[idx])
                    time.sleep(2)
                control.inSignal_1(self)
            control.inSignal_2(self)
        control.inSignal_term(self)
                # Signal Serpent (SIGUSR2) #
#            else:
#                # Signal Serpent (SIGUSR2) #
#                control.inSignal_2(self)

            # Signal Serpent (SIGTERM) #
#            if not self.simulating:
#                break
        

 #           control.outSignal_term(self)
        

    # Wait for signal 
    def waitSignal(self):
            self.sleeping = True
            while self.sleeping:
                # Sleep for two seconds
                time.sleep(2)
                # Open file to check if we got a signal
                fin = open(self.outSignal,'r')
                # Read line 
                line = fin.readline()
                # Close file 
                fin.close()
                # Check signal
                if int(line) != -1:
                 #   print (line)
                    if int(line) == signal.SIGUSR1:
                        # Got the signal to resume
                        self.sleeping = False
        
                    elif int(line) == signal.SIGUSR2:
                        # Got the signal to move to next time point
                        self.iterating = False
                        self.sleeping = False

                    elif int(line) == signal.SIGTERM:
                        # Got the signal to end the calculation
                        self.iterating = False
                        self.sleeping = False
                        self.simulating = False
                    else:
                        # Unknown signal
                        print ("\nUnknown signal read from file, exiting\n")
                        # Exit
                        quit()                
        
                    # Reset the signal in the file
                    file_out = open(self.outSignal,'w')
                    file_out.write('-1')
                    file_out.close()
 
    # Read power distribution #
    def pRead(self,outfile_step):
        outfile = 'OUT_'+outfile_step
        # Reset power distribution
        self.Power = []
        self.p_Error = []
        # Open output file
        file_in = open(outfile,'r')
        # Loop over output file to read power distribution
        for line in file_in:
                if line.split()[0] == 1:
                    # Split line to values
                    value = line.split()
                    # Store power 
                    self.Power.append(float(value[-2]))
                    self.p_Error.append(float(value[-1]))
        file_in.close()
    # check convergency
    def tConverge(self,epsi):
        convT = []
        with open (self.conv_T,'r') as f:
            for line in f:
                if 'T_eps(idx)' in line:
                    data = line.split()[-1]
                    # get rid of ;
                    data = data.rstrip(';')
                    convT.append(float(data))
        f.close()
        if convT[-1] < epsi:
            print ('Temperature Converged')
            conv = True
        else:
            conv = False

        return conv
    # send signal to com.in
    def inSignal_1(self):   
        file_out = open(self.inSignal,'w')
        file_out.write(str(int(signal.SIGUSR1)))
        file_out.close()
    def inSignal_2(self):   
        file_out = open(self.inSignal,'w')
        file_out.write(str(int(signal.SIGUSR2)))
        file_out.close()
    def inSignal_term(self):   
        file_out = open(self.inSignal,'w')
        file_out.write(str(int(signal.SIGTERM)))
        file_out.close()
    # send signal to com.out
    def outSignal_1(self):   
        file_out = open(self.outSignal,'w')
        file_out.write(str(int(signal.SIGUSR1)))
        file_out.close()
    def outSignal_2(self):   
        file_out = open(self.outSignal,'w')
        file_out.write(str(int(signal.SIGUSR2)))
        file_out.close()
    def outSignal_term(self):   
        file_out = open(self.outSignal,'w')
        file_out.write(str(int(signal.SIGTERM)))
        file_out.close()
        self.interating = True
        self.sleeping = True


    # write generated ifc file to serpent read ifc file
    def writeIFC(ele_fArray):
        ifcPre = ele_fArray.split('_')[0]
        ifcFile = ifcPre+'_'+'file'
        ifc_curr = ele_fArray

        with open (ifc_curr, 'r') as f:
            inContent = f.readlines()
        f.close()

        with open (ifcFile, 'w') as f:
            for line in inContent:
                f.write(line)
        f.close()
    # add name of convergence file to class
    def convFile(self):
        ifcPre = self.fArray[0].split('_')[0]
        ifcFile = ifcPre+'_'+'file'
        # add ifcFile to data
        self.ifcFIle = self.ifcFile + ifcFile
        self.conv_P = self.conv_P + ifcFile + '_Pconv0.m'
        self.conv_T = self.conv_T + ifcFile + '_Tconv0.m'

    # generate array of FEMAXI time step
    def fTimeGen(self):
        time = []
        for name in self.fArray:
            timeStr = name.split('_')[-1]
            timePoint= float(timeStr.replace('d','.'))
            time.append(timePoint)
        self.fTime = self.fTime + time

    # Python code to remove duplicate elements 
    def reDup(self): 
        time = [] 
        for ele in self.fTime: 
            if ele not in time: 
                time.append(ele) 
        self.fTime = time

    # generate array of Serpent time step
    def sTimeGen(self):
        with open (self.sInput,'r') as f:
            content = f.readlines()
        f.close()
        # get timestep keyword
        for line in content:
            if 'nps' in line and not line.startswith('%'):
                timeKey = line.split()[-1]
        for line in content:
            if 'tme' in line and timeKey in line and not line.startswith('%'):
               # print (line)
                tType = int(line.split()[2])
                if tType != 1: 
                    nBins = int(line.split()[3])
                    sTime = float(line.split()[4])
                    eTime = float(line.split()[5])
        if tType == 2:
            step =  (eTime - sTime)/nBins
#            sTime = sTime + step
            eTime = eTime + step
            time = np.arange(sTime,eTime,step)

        time = list(time)
        self.sTime = self.sTime + time
       # self.sTime = self.sTime + sTime
#name = 'tf3'
#ifctype = 5
#axial = [0.0,400.0]
#angular = [1,0.0,360.0]
#radial = [0.0]
#segPin = 1
#segName = ['pina','pinb','pinc']
##pin1 = pin(name,ifctype,axial,angular,radial)
#pin1 = pinData(name,ifctype,segPin)
#pin1.addAxial(axial)
#pin1.addNcrack(angular[0])
#nCrack = pinData.get_n_crack(pin1)
#
#bool_segPin = pinData.segPin_bool(pin1)
#if pin1.segPin_bool():
#    pin1.addSeg(segName)
#
#infile = 'inputData'
#
#inData = inReader(infile)
#inData.addName()
#inData.addIFC()
#inData.addSegArray()
#inData.addAxial()
#inData.addAngular()
#inData.addRadial()
