#Elliott Mueller and Fenfang Wu, GPS Caltech, 2022
import numpy as np
from numpy import *
from tqdm import tqdm
from scipy import optimize
import pylab as py
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.integrate import odeint
ratio_VPDB=0.01118 #reference ratio
ratio_naturalC = ratio_VPDB
abun_C13_natural=ratio_naturalC/(ratio_naturalC+1)  # fractional abundance of rare isotope (here 13C, but can be any isotope) 
abun_C12_natural=1-abun_C13_natural

def GetProp(conc):  #INPUTS: isotopologue abundances for a given molecule OUTPUTS: Array of labelling proportions M, M+1, M+2, etc. where M is the molecular weight. Size = (1 x number of atoms + 1)
    nm=len(conc)
    nC=int(np.log2(nm))
    total=sum(conc)
    M=np.zeros(nC+1)
    for im in range(nm):
        nC13=0
        n=im
        while(n):
            nC13 += n&1
            n >>=1
        M[nC13]+=conc[im]
    
    return M/total
    
def siteSpecific(sname, nC,conc): #INPUTS: sname = name of molecule, nC = number of atoms in molecule, conc = isotopologue abundances for molecule   OUTPUTS: Array containing [name of compound, site specific isotope composition in delta notation relative to the standard ratio defined above.]
    output = [sname]
    for iC in range(nC):
        cumu_ratio=GetSingleRatio(conc,iC)
        cumu_dC=GetDelta(cumu_ratio)
        output.append("C{}: {} ".format(iC+1, cumu_dC))
    return output
    #print(GetProp(conc))

def GetSingleRatio(conc, iC):#INPUTS: conc = isotopologue abundances for molecule, iC is the atomic number of interest OUTPUTS: Single isotope ratio value for the atomic site (iC + 1). If iC = 1, the isotope ratio of atomic site 2 will be reported.
    nm=len(conc)
    nC=int(np.log2(nm))
    ib=nC-(iC+1)  #the bit number, 5,4,3,2,1,0 from left to right
    total_C13=0
    total_C12=0
    for im in range(nm):
        if (im>>ib)&1==1:
            total_C13+=conc[im]
        else:
            total_C12+=conc[im]
    return total_C13/total_C12
    

def GetDelta(ratio): # INPUTS: ratio = isotope ratio OUTPUTS: delta value against reference ratio above ("ratio_VPDB")
    delta=round((ratio/ratio_VPDB-1)*1000*1000)/1000.
    return delta

def Intramolecular(nC,intra_structure,initial):#INPUTS: nC = number of carbons in molecule, intra_structure = user defined vector from atomic site 1 to nC with delta values (‰) for each site in the molecule, initial = initial concentration of the molecule. OUTPUTS: conc_t0 = starting isotopologue abundance array based on initial intramolecular composition (size = 1 x 2^nC)
    nm=2**nC
    conc_t0=np.repeat(initial,nm)
    n_C13=np.zeros(nm)
    n_C12=np.zeros(nm)
    abun_C13=0
    abun_C12=0
    abun_C13_natural=ratio_naturalC/(ratio_naturalC+1)
    intra_structure_ratio = (intra_structure/1000+1)*0.01118
    intra_structure_abun = intra_structure_ratio/(intra_structure_ratio+1)
    
    for im in range(nm):
        for ib in range(nC-1,-1,-1):#left to right, high bit to low bit
            abun_C13 = intra_structure_abun[(nC-1)-ib]
            abun_C12 =1-abun_C13
            if (im>>ib)&1==1:
                conc_t0[im] *= abun_C13
                n_C13[im]+=1
            else:
                conc_t0[im]*=abun_C12
                n_C12[im]+=1
                
    return conc_t0


def GetFullRates(nC, singlerates): #INPUTS: nC = number of atoms in molecule, singlerates = an array with monoisotopic and single substituted rate constants, defined in the NB file. OUTPUTS: array of reaction rate constants for every isotopologue in a reactant which represents the user-defined fractionation factors.
    if len(singlerates)==1: #no fractionation
        return np.repeat(singlerates,2**nC) # if singlerates (defined
    else:
        nm = 2**nC
        full_rates=np.repeat(1.0, nm)
        full_rates[0]=singlerates[0]
        for im in range(nm):
            nC13=0
            for iC in range(nC-1,-1,-1):
                if im>>iC &1 == 1: #C13 at iC position
                    #breakpoint()
                    full_rates[im]*=singlerates[nC-iC]
                    nC13+=1
                else:
                    full_rates[im]*=1
            if im>0:
                full_rates[im]/=singlerates[0]**(nC13-1)
        return full_rates

def synthesis(conc_a,conc_b,rate_abx): # A + B -> C, INPUTS: conc_a = isotopologue abundances of Reactant A, conc_b = isotologue abundances of Reactant B, rate_abx = full rates array (size = 1 x 2^(nC_A + nC_B)) OUTPUTS: created_x = isotopologue vector of product,
    ma = len(conc_a)
    iC=int(np.log2(ma))
    mb = len(conc_b)
    conc_a_expanded = np.repeat(conc_a,mb) #expand concentration array of Reactant A
    conc_b_expanded = np.tile(conc_b,ma) #expand concentration array of Reactant A
    created_x = conc_a_expanded*conc_b_expanded*rate_abx #run reaction to calculate amount of product made
    reacted_a, reacted_b=breakdown(created_x,iC) #calculate amount of reactant A and B used in the reaction by running the breakdown function on the product

#    if np.any(np.isinf(created_x)):
 #       breakpoint()
    return (created_x, reacted_a, reacted_b)


def breakdown(reacted_x, nC_a): #INPUTS: array of isotopologues for a given reactant, nC_a = number of atoms in reactant A. Tells QIRN where to break the reactant in two. OUTPUTS: Calculated isotopologue arrays for creation of two products (created_a, created b)
    mx =len(reacted_x)
    ma=2**nC_a
    mb=round(mx/ma)
    reacted_x.shape = (ma, mb)
    created_a=reacted_x.sum(1) #sum each rows, in total ma rows
    created_b=reacted_x.sum(0) #sum each column, in total mb columns
    reacted_x.shape =(ma*mb)
    return(created_a, created_b)

def MapCchain_index(old_ibit, new_ibit): #This function creates a vector that tells QIRN how to rearrange the atomic positions of a molecule while retaining mass balance of 13C and 12C. INPUTS: old_ibit = reference barcode for the molecule [nC -1 -> 0], new_ibit = barcode telling QIRN how to rearrange the atomic positions. OUTPUTS: new_index = an array of integers where each integer represents the index of the array of isotopologues representing the transformed molecule.
    #nm=len(conc)
    nC=len(old_ibit)
    nm=2**nC;
    new_index=zeros(nm,dtype='int')
    for im in range(nm):
        new_im=0
        for iC in range (nC):
            new_im+=((im>>old_ibit[iC]) & 1)<<new_ibit[iC]
        new_index[new_im]=im
        #new_conc[new_im]=conc[im]
        #print (im, new_im)
    #To use:    
    #new_conc=conc[new_index]
    return new_index

def GetTotalRatio(nC, conc_x): #get a compound-specific isotope ratio INPUT: nC = number of atoms in the molecule, conc_x = isotopologue array for that molecule
    nm=2**nC
    n_C13=zeros(nm)
    n_C12=zeros(nm)
    total_C12=0
    total_C13=0
    for im in range(nm): 
        for iC in range(nC-1, -1,-1):
            if im>>iC &1==1:
                n_C13[im]+=1
            else:
                n_C12[im]+=1
        total_C12 += conc_x[im]*n_C12[im]
        total_C13 += conc_x[im]*n_C13[im]
    return total_C13/total_C12


if __name__=="__main__":
    print('util')
    
def QIRN(intermediates,networkbuilder,reactiondatabase, time, dt): #Main function for running the numerical QIRN model
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import csv
    from numpy import zeros,array,log,concatenate
    from copy import deepcopy


    intmed={} #a dictionary
    metanetwork = {} #a dictionary that collates all the relevant information from the NB, RD and IntMedDatabase to construct the user-defined network
    substrates = [] #list of substrates used in the network
    reactions = [] #list of reactions based on their reaction names from RD
    reservoirs = [] #list of molecules that are treated as reservoirs and do not change concentration or isotope composition in time
    EClist = [] #list of enzymes based on their Reaction ID's from the RD and NB

#SETTING INITIAL CONDITIONS
    with open(intermediates,newline='') as imfile: #open IntMedDatabase
        imreader = csv.reader(imfile)
        next(imreader)  #get rid of the first row
        for row in imreader: #populate isotopologue arrays for the initial concentrations and isotope compositions of the intermediates.
            name=row[0].strip()
            initial=float(row[2])
            reservoir = row[4]
            nC = int(row[1])
            if initial>0:
                if row[3] == '':
                    structure = np.array(np.zeros(nC),dtype = 'f')
                else:
                    structure =[float(i) for i in row[3].split(',')]
                    structure = np.array(structure, dtype = 'f')
                conc = Intramolecular(nC,structure,initial)
                conc_init = Intramolecular(nC,structure,initial)
                
            else:
                conc=np.zeros(2**nC)
                conc_init = np.zeros(2**nC)
            intmed[name]=[nC,initial,conc,np.zeros(2**nC),conc_init,reservoir] #intmed is a dictionary that keeps track of changes in concentration of every isotopologue in the network. Each row is a molecule. nC = number of atoms, initial = initial concentration set by in the IntMedDB, conc = integrate isotopologue abundance array, fourth index = instantaneous change of the molecule after all the reactions of a timestep have run (this can be + or -), conc_init = the initial isotopologue abundance array, reservoir = if on, put this molecule in the reservoir list.
            if intmed[name][4][0] >0 and intmed[name][5]!='': #create list of reservoir substrates
                reservoirs.append(name)


#POPULATING METANETWORK AND CONSTRUCTING USER-DEFINED NETWORK
    with open(networkbuilder, newline='') as networkfile: #open NB
            networkreader = csv.reader(networkfile)
            next(networkreader)

            for row in networkreader: #Populate metanetwork dictionary using Reaction ID (EC number, etc.) as the key. Loop through NB to do this for every reaction in the in network.
                kf = row[1] #forward reaction rate
                kr = row[2] #reverse reaction rate
                EC = row[0].strip()
                metanetwork[EC] = [kf, kr]
                metanetwork[EC][2:10] = np.zeros((7,),dtype=int)
                metanetwork[EC][5:6] = [kf,kr]
                EClist.append(EC) #create list of reactions based on Reaction ID


    with open(reactiondatabase, newline='') as reactiondatabasefile: #Open RD and populate metanetwork rows with the reactant and product names for each reaction, as well as any necessary MapChain arrays to rearrange molecules during reactions.
        databasereader = csv.reader(reactiondatabasefile)
        next(databasereader)
        for row in databasereader:
            try:
                metanetwork[row[0]]
            except:
                continue

            EC = row[0]
            reactantA = row[2]
            reactantB = row[3]
            productC = row[4]
            productD = row[5]
            reaction = row[1]

            if len(row[8])>0: #test whether a reaction requires molecular rearrangemennt
                old_ibit=[int(i) for i in row[8].split(',')] #get old and new ibits from the RD, these are always the same for a given reaction.
                newf_ibit=[int(i) for i in row[6].split(',')]
                newr_ibit=[int(i) for i in row[7].split(',')]
                findex = MapCchain_index(old_ibit,newf_ibit) #create rearrangement array for the forward reaction
                rindex = MapCchain_index(old_ibit,newr_ibit) #create rearrangement array for the reverse reaction
                metanetwork[EC][2:5] = [reactantA, reactantB, productC, productD] #populate metanetwork with Reactant A, B and Product C, D (strings)
                metanetwork[EC][8] = findex
                metanetwork[EC][9] = rindex
                metanetwork[EC][10] = reaction
                #print(reaction)
            else:
                metanetwork[EC][2:5] = [reactantA, reactantB, productC, productD]
                metanetwork[EC][10] = reaction
                #print(reaction)
  

##CONVERT MASS ACTION RATE LAWS TO VECTORS WITH ASSOCIATED ISOTOPE EFFECTS APPLIED
    
    sub = []
    with open(intermediates, newline='') as imfile:
            imreader = csv.reader(imfile)
            next(imreader)  #get rid of the first row
            for row in imreader:
                sub.append(row[0])
    
    substrates = []
    for i in EClist:
        for j in sub:
            if j in metanetwork[i][2:6] and j not in substrates:
                substrates.append(j)
                
    with open(networkbuilder, newline='') as networkfile:
        networkreader = csv.reader(networkfile) 
        next(networkreader)
        for row in networkreader:
            try:
                x = metanetwork[row[0]][10]
            except IndexError:#check whether reactions in NB match reactions in the RD
                print('Looks like this reaction does not match anything in the database: ', row[0])
                #sys.exit()
            EC = row[0]
            ir = metanetwork[EC]
            reactions.append(ir[10])
            k_for = ir[0]
            k_rev = ir[1]
            nC=intmed[ir[2]][0]
            if len(ir[3])>0:
                nC+=intmed[ir[3]][0]
            if len(row[3]) > 0: #convert site specific alpha values for the forward reaction from the NB into an array of rate constants (size = 1 x 2^N+1 where N is the number of carbon positions)
                k_for = np.repeat(ir[0], nC+1)
                k_for_Csites=[int(i) for i in row[3].split(',')]
                k_for_KIEs=[float(i) for i in row[4].split(',')]
                k_for_weights = np.ones(nC)

                j = 0
                for i in k_for_Csites:
                    k_for_weights[i-1] = k_for_KIEs[j]

                    j = j +1
                k_for_weights = np.array(k_for_weights, dtype = 'f')
                k_for = np.array(k_for, dtype = 'f')
                k_for[1::] = k_for[1::] * k_for_weights

            else:
                k_for=[float(i) for i in ir[0].split(',')]

            if len(row[5]) > 0: #convert site specific alpha values for the reverse reaction from the NB into an array of rate constants (size = 1 x 2^N+1 where N is the number of carbon positions)
                k_rev = np.repeat(ir[1], nC+1)
                k_rev_Csites=[int(i) for i in row[5].split(',')]
                k_rev_KIEs=[float(i) for i in row[6].split(',')]
                k_rev_weights = np.ones(nC)    

                j = 0 
                for i in k_rev_Csites:
                    k_rev_weights[i-1] = k_rev_KIEs[j]
                    j = j + 1    
                k_rev_weights = np.array(k_rev_weights, dtype = 'f')
                k_rev = np.array(k_rev, dtype = 'f')
                k_rev[1::] = k_rev[1::] * k_rev_weights
            else:
                k_rev=[float(i) for i in ir[1].split(',')]

            metanetwork[EC][0]=GetFullRates(nC,k_for)*dt #populate
            metanetwork[EC][1]=GetFullRates(nC,k_rev)*dt

    isocitr_channel=2  # Symmetry toggle for the actonitase reaction (citrate -> isocitrate) 1 = C4 channel receives hydroxyl group, 2 = Purdue's C2 site receives hydroxyl group, 0 = both C2 and C4 receive hydroxyl group (split in half)
    #breakpoint()


    t =time #total reaction time in second
    skip = 50 #cadence for bookkeeping of the reaction fluxes and concentrations. Raising this value will make the time steps run faster but lose resolution on the changes in flux, concentration and isotope composition with time. This does NOT impact the number of time steps run. Only the size of the dictionaries that keep track of the network's parameters in time.
    MProportions = np.zeros((12,len(substrates)))
    t_steps=int(t/dt) #number of timesteps to take
    
    concentration_tracker = np.zeros((len(substrates),int(t_steps/skip)+1)) #an array which holds the total concentration for every molecule (sum of all isotopologues) in the NB through time. Time points only taken at a certain cadence set by the variable 'skip'
    isotope_tracker = np.zeros((len(substrates),int(t_steps/skip)+1)) #an array which holds the compound-specific isotope composition for every molecule in the NB through time. Time points only taken at a certain cadence set by the variable 'skip'
    flux_tracker = np.zeros((len(reactions),int(t_steps/skip)+1)) #an array which holds the net (forward - reverse) fluxes for every reaction in the NB through time. Time points only taken at a certain cadence set by the variable 'skip'
    
    with open(networkbuilder, newline='') as networkfile: #create list of reactions from reaction ID's because this was lost in previous portion of the code.
        networkreader = csv.reader(networkfile)
        next(networkreader)
        i = 0
        EC_list = {}
        for row in networkreader:
            EC_list[i] = row[0]
            i = i + 1

##WHERE BOX MODEL BEGINS
    for it in tqdm(range(t_steps+1)):
        it_check = it/skip #checks whether it is a time step for which it needs to record parameters
        m = 0
        for row in EC_list: #loops through each row of the metanetwork (represents a reaction)
            ir = metanetwork[EC_list[row]]
        # Reaction A -> C
            if ir[3]==''and ir[5]=='': #Is it a simple transformation?
                reactedA=intmed[ir[2]][2]*ir[0] #forward
                createdC=reactedA
                reactedC=intmed[ir[4]][2]*ir[1]#reverse
                createdA=reactedC
                
                if ir[8].size>1 and ir[9].size>1:  #does this need to be remapped?
                    createdC = reactedA[ir[8]]
                    createdA = reactedC[ir[9]]
                if ir[2]=='citrate': #citrate -> isocitrate symmetry
                    if isocitr_channel==0: #two positions to add OH
                        createdC4 = 0.5*reactedA #isocitrate C4-OH
                        createdC2 = createdC4[ir[8]]
                        createdC=createdC2+createdC4
                    elif isocitr_channel==2: #isocitrate C2-OH
                        createdC = reactedA[ir[8]]
                if ir[2]=='fumarate': #fumarate->malate symmetry
                    createdC2 = 0.5*reactedA
                    createdC3 = createdC2[ir[8]]
                    createdC=createdC2+createdC3

                intmed[ir[2]][3]+=createdA - reactedA #calculate change in concentration for Reactant A and Product C and add it to the instantaneous change index in the 'intmed' dictionary
                intmed[ir[4]][3]+=createdC - reactedC
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA) #flux calculation based on the change in the total concentration of reactant A (sum of isotopologues).


        # Breakdown reaction  A->C+D
            elif ir[3]=='' and ir[5]!='': #is it a breakdown function?
                reactedA=intmed[ir[2]][2]*ir[0] #forward
                nC_C = intmed[ir[4]][0]

                createdC,createdD=breakdown(reactedA,nC_C)
                createdA,reactedC,reactedD=synthesis(intmed[ir[4]][2],intmed[ir[5]][2],ir[1])#reverse
                
                
                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    reactedM = reactedA[ir[8]]
                    createdC,createdD=breakdown(reactedM,nC_C)
                    createdA = createdA[ir[9]]

                intmed[ir[2]][3] +=createdA - reactedA #calculate change in concentration for Reactant A Product C and Product D and add it to the instantaneous change index in the 'intmed' dictionary
                intmed[ir[4]][3] +=createdC - reactedC
                intmed[ir[5]][3] +=createdD - reactedD
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)  #flux calculation based on the change in the total concentration of reactant A (sum of isotopologues).



        # Synthesis reaction  A+B->C
            elif ir[3]!='' and ir[5]=='':
                createdC,reactedA,reactedB=synthesis(intmed[ir[2]][2],intmed[ir[3]][2],ir[0])#forward
                
                reactedC = intmed[ir[4]][2]*ir[1]
                createdA,createdB=breakdown(reactedC,intmed[ir[2]][0]) #reverse

                if ir[4]=='citrate': #acecoa+oxaace->citrate
                    createdC = createdC[ir[8]]
                    reactedC = reactedC[ir[9]]

                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    createdC = createdC[ir[8]]
                    reactedM = reactedC[ir[9]]
                    createdA,createdB=breakdown(reactedM,intmed[ir[2]][0])

                intmed[ir[2]][3] +=createdA -reactedA #calculate change in concentration for Reactant A Reactant B and Product C and add it to the instantaneous change index in the 'intmed'
                intmed[ir[3]][3] +=createdB - reactedB
                intmed[ir[4]][3] +=createdC - reactedC
                if it_check.is_integer(): #check if it is a time step where you need to record.
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA) #flux calculation based on the change in the total concentration of reactant A (sum of isotopologues).
                
                    if ir[2] == ir[3]: #check the stoichiometry of the reaction. If two of the same molecule are condensed, the flux is twice as high for that molecule but the product creation flux is half that. (only the reaction consumption flux is reported in the flux_tracker array)
                        flux_tracker[m][int(it_check)] *=2

        # A + B --> C + D
            elif ir[2]!='' and ir[5]!='':
                createdM, reactedA,reactedB = synthesis(intmed[ir[2]][2],intmed[ir[3]][2],ir[0]) #make intermediate M by synthesis function on A and B
                rev_createdM,reactedC,reactedD = synthesis(intmed[ir[4]][2],intmed[ir[5]][2],ir[1]) #make intermediate M for the reverse reaction by synthesis function on C and D
                #breakpoint()
                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    createdM = createdM[ir[8]] #remap the intermediate using the index transformation array in the metanetwork
                    rev_createdM = rev_createdM[ir[9]] #remap the intermediate using the index transformation array in the metanetwork

                elif ir[5] == 'co2_als': #special case for the reaction pyruvate + pyruvate -> acetolactate + CO2 (april 2022, EPM)
                    createdM = createdM[pyr_pyr_index]*0.5
                    rev_createdM = rev_createdM[pyr_pyr_index]

                nC_C=intmed[ir[4]][0] #figure out how big the Product C and Product D are for the forward reaction and break the intermediate at this atomic site to create the final products.
                createdC,createdD=breakdown(createdM,nC_C)
                nC_A=intmed[ir[2]][0]
                createdA,createdB = breakdown(rev_createdM,nC_A) #figure out how big the Reactant A and Reactant B are for the reverse reaction and break the intermediate at this atomic site to create the final products.

                intmed[ir[2]][3] +=createdA - reactedA #calculate change in concentration for Reactant A Reactant B, Product C, and Product D and add it to the instantaneous change index in the 'intmed'
                intmed[ir[3]][3] +=createdB - reactedB
                intmed[ir[4]][3] +=createdC - reactedC
                intmed[ir[5]][3] +=createdD - reactedD
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)
                
                    if ir[2] == ir[3]: #check the stoichiometry of the reaction. If two of the same molecule are condensed, the flux is twice as high for that molecule but the product creation flux is half that. (only the reaction consumption flux is reported in the flux_tracker array)
                        flux_tracker[m][int(it_check)] *=2
                    
            m = m+1
    #changing concentration of intermediates
        i = 0
        
        if sum(intmed[row][2]) > 0 and it_check.is_integer(): #check whether substrates have a non-zero concentration and whether the time step is being recorded
            for row in substrates:
                concentration_tracker[i][int(it_check)]=sum(intmed[row][2]) #track concentration (sum of isotopologues)
                isotope_track = GetTotalRatio(intmed[row][0],intmed[row][2]) #calculate compound-specific isotope ratio for each substrate
                isotope_tracker[i][int(it_check)] = GetDelta(isotope_track) #track compound-specific delta value for each substrate
        
        for row in substrates: #at every timestep update integrated concentration and reset the instantaneous change in concentration to zero. Set reservoir molecules to their initial concentration and composition
            intmed[row][2] += intmed[row][3]
            intmed[row][3] = np.zeros(2**intmed[row][0])
            if row in reservoirs:
                intmed[row][2][:] = intmed[row][4][:]
            
            i = i+1
        
    #const supply of the following reagents
            


        if it+2==t_steps:
            prev_intmed = deepcopy(intmed)
        
        flux_tracked = {}
            
    for i in range(len(reactions)):
        flux_tracked[reactions[i]] = flux_tracker[i][:] #create array of fluxes to be used in the GUI for data visualization
    
    ## NETWORK DESCRIPTIONS:
    print('Number of reactions: ',len(EClist))     # Number of reactions
    print('Number of substrates: ',len(substrates)) # Number of substrates
    isotopologues = 0
    
    for row in substrates:
        isotopologues += 2**(intmed[row][0]) #Number of isotopologues
    print('Number of isotopologues: ', isotopologues)

    return intmed, substrates, concentration_tracker, reactions, flux_tracker, metanetwork, EClist,flux_tracked, reservoirs, isotope_tracker, skip #return values to the GUI for data visualization

def emprint(compound,intmed): #prints a compound's name followed by its site specific composition
    value = siteSpecific(compound, intmed[compound][0], intmed[compound][2])
    return value
