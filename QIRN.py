#Fenfang Wu, GPS Caltech, 2020
#!env python3
import numpy as np
from numpy import *
from tqdm import tqdm
from scipy import optimize
import pylab as py
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#ratio_naturalC = 0.01121437
ratio_VPDB=0.01118  #C13/C12
ratio_naturalC = ratio_VPDB
abun_C13_natural=ratio_naturalC/(ratio_naturalC+1)  # C13% over total
abun_C12_natural=1-abun_C13_natural

#delta_labeled = 10000. #per mil
#ratio_labeled = (delta_labeled/1000+1)*ratio_VPDB
ratio_labeled= 99.99 #C13/C12
abun_C13_labeled=ratio_labeled/(ratio_labeled+1)
abun_C12_labeled =1-abun_C13_labeled

#GITHUB CHECK
def GetProp(conc):
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
    
    

def write2file(ofile, sname, created, conc_t):
    nC=int(round(np.log2(len(created))))

    ofile.write("\n{0}: \n Inst amount:  {1} moles  \n".format(sname,sum(created)))          
    #for iC in range(nC):
        #ofile.write("C{}:".format(iC+1))
        #inst_ratio=GetSingleRatio(created,iC)
        #cumu_ratio=GetSingleRatio(conc_t,iC)
        #inst_dC=GetDelta(inst_ratio)
        #cumu_dC=GetDelta(GetSingleRatio(conc_t,iC))
        #ofile.write("{0:.12f}   {1:.12f}   {2:.3f}     {3:.3f}\n".format(inst_ratio, cumu_ratio, inst_dC,cumu_dC))


def siteSpecific(sname, nC,conc):
    #print(sname)
    output = [sname]
    for iC in range(nC):
        cumu_ratio=GetSingleRatio(conc,iC)
        cumu_dC=GetDelta(cumu_ratio)
        #print("C{}: {} ".format(iC+1, cumu_dC))
        output.append("C{}: {} ".format(iC+1, cumu_dC))
    return output
    #print(GetProp(conc))



       
def GetSingleRatio(conc, iC):#iC is the labeled C (0,1,2,3,4,5 from left to right)
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
    

def GetDelta(ratio): # C delta value
    delta=round((ratio/ratio_VPDB-1)*1000*1000)/1000.
    return delta

def Intramolecular(nC,intra_structure,initial):#Get conc for each isotopomer, iC is the labeled C... intrastructure is a vector [C1->C6] in delta permil
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

def GetConc(iC,nC):#Get conc for each isotopomer, iC is the labeled C
    nm=2**nC
    conc_t0=np.repeat(1.,nm)
    n_C13=np.zeros(nm)
    n_C12=np.zeros(nm)
    abun_C13=0
    abun_C12=0
    for im in range(nm):
        for ib in range(nC-1,-1,-1):#left to right, high bit to low bit
            if ib==(nC-iC-1):
                abun_C13 = abun_C13_labeled
            else:
                abun_C13 = abun_C13_natural
            abun_C12 =1-abun_C13
                       
            if (im>>ib)&1==1:
                conc_t0[im] *= abun_C13
                n_C13[im]+=1
            else:
                conc_t0[im]*=abun_C12
                n_C12[im]+=1
                
           
    #print(sum(conc_t0))
    return conc_t0

def GetFullRates(nC, singlerates): #single rate is an array with mobsnoisotopic and single suituted rate constants, nC+1 numbers
    if len(singlerates)==1: #no fractionation
        return np.repeat(singlerates,2**nC)
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

def synthesis(conc_a,conc_b,rate_abx): # A + B -> X
    ma = len(conc_a)
    iC=int(np.log2(ma))
    mb = len(conc_b)
    conc_a_expanded = np.repeat(conc_a,mb)
    conc_b_expanded = np.tile(conc_b,ma)#no need to be normalized
    created_x = conc_a_expanded*conc_b_expanded*rate_abx
    reacted_a, reacted_b=breakdown(created_x,iC)

    if np.any(np.isinf(created_x)):
        breakpoint()
    return (created_x, reacted_a, reacted_b)

#arr_x is the reacted amount of x
def breakdown(reacted_x, nC_a): #break the array x to array a and b
    mx =len(reacted_x)
    ma=2**nC_a
    mb=round(mx/ma)
    reacted_x.shape = (ma, mb)
    created_a=reacted_x.sum(1) #sum each rows, in total ma rows
    created_b=reacted_x.sum(0) #sum each column, in total mb columns
    reacted_x.shape =(ma*mb)
    return(created_a, created_b)

def MapCchain(conc, old_ibit, new_ibit):
    nm=len(conc)
    nC=len(old_ibit)
    new_conc=zeros(nm)
    for im in range(nm):
        new_im=0
        for iC in range (nC):
            new_im+=((im>>old_ibit[iC]) & 1)<<new_ibit[iC]
        new_conc[new_im]=conc[im]
        #print (im, new_im)   
    return new_conc

def MapCchain_index(old_ibit, new_ibit):
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

def GetTotalRatio(nC, conc_x): #get a compound's total_C13/total_C12
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
    
def QIRN(intermediates,networkbuilder,reactiondatabase, time, dt):
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import csv
    from numpy import zeros,array,log,concatenate
    from copy import deepcopy


    intmed={} #a dictionary
    metanetwork = {} #a dictionary
    substrates = []
    flux_tracker ={}
    reactions = []
    reservoirs = []
    EClist = []


    with open(intermediates,newline='') as imfile:
        imreader = csv.reader(imfile)
        next(imreader)  #get rid of the first row
        for row in imreader:
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
                
                #conc=GetConc(-1,nC)
                
                #if name=='3-phosphoglycerate':
                   # conc=GetConc(-1,nC)
                   # conc_3PGA[:] = conc
                #if name =='glyoxylate':
                 #   conc=GetConc(1,nC)
                  #  conc_glyoxylate[:] = conc
            else:
                conc=np.zeros(2**nC)
                conc_init = np.zeros(2**nC)
            intmed[name]=[nC,initial,conc,np.zeros(2**nC),conc_init,reservoir]
            if intmed[name][4][0] >0 and intmed[name][5]!='':
                reservoirs.append(name)



    with open(networkbuilder, newline='') as networkfile:
            networkreader = csv.reader(networkfile)
            next(networkreader)

            for row in networkreader:
                kf = row[1]
                kr = row[2]
                EC = row[0].strip()
                metanetwork[EC] = [kf, kr]
                metanetwork[EC][2:10] = np.zeros((7,),dtype=int)
                metanetwork[EC][5:6] = [kf,kr]
                EClist.append(EC)


    with open(reactiondatabase, newline='') as reactiondatabasefile:
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

            if len(row[8])>0:
                old_ibit=[int(i) for i in row[8].split(',')]
                newf_ibit=[int(i) for i in row[6].split(',')]
                newr_ibit=[int(i) for i in row[7].split(',')]
                findex = MapCchain_index(old_ibit,newf_ibit)
                rindex = MapCchain_index(old_ibit,newr_ibit)
                metanetwork[EC][2:5] = [reactantA, reactantB, productC, productD]
                metanetwork[EC][8] = findex
                metanetwork[EC][9] = rindex
                metanetwork[EC][10] = reaction
                #print(reaction)
            else:
                metanetwork[EC][2:5] = [reactantA, reactantB, productC, productD]
                metanetwork[EC][10] = reaction
                #print(reaction)
  

 #   e = 0
  #  for row in EClist:
   #     metanetwork[row][0] = str(optval_f[e])
    #    metanetwork[row][1] = str(optval_r[e])
     #   e +=1
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
            except IndexError:
                print('Looks like this reaction does not match anything in the database: ', row[0])
                sys.exit()
            EC = row[0]
            ir = metanetwork[EC]
            reactions.append(ir[10])
            k_for = ir[0]
            k_rev = ir[1]
            nC=intmed[ir[2]][0]
            if len(ir[3])>0:
                nC+=intmed[ir[3]][0]
            if len(row[3]) > 0:
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

            if len(row[5]) > 0:
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

            #if len(row[3])==0:
               # k_for=[float(i) for i in ir[0].split(',')]
               # k_rev=[float(i) for i in ir[1].split(',')]
                
            metanetwork[EC][0]=GetFullRates(nC,k_for)*dt
            metanetwork[EC][1]=GetFullRates(nC,k_rev)*dt

    PK_equilibration = 0
    isocitr_channel=2  #1 means Alex's C4 channel, 2 means Purdue's C2 channel, 0 means both channel
    #breakpoint()


    t =time #total reaction time in second
    skip = 50
    MProportions = np.zeros((12,len(substrates)))
    t_steps=int(t/dt)
    
    concentration_tracker = np.zeros((len(substrates),int(t_steps/skip)+1))
    isotope_tracker = np.zeros((len(substrates),int(t_steps/skip)+1))
    flux_tracker = np.zeros((len(reactions),int(t_steps/skip)+1))
    
    with open(networkbuilder, newline='') as networkfile:
        networkreader = csv.reader(networkfile)
        next(networkreader)
        i = 0
        EC_list = {}
        for row in networkreader:
            EC_list[i] = row[0]
            i = i + 1

        
    ##WHERE BOX MODEL BEGINS
    for it in tqdm(range(t_steps+1)):
        it_check = it/skip
        m = 0
        for row in EC_list:
            ir = metanetwork[EC_list[row]]
        # Reaction A -> C
            if ir[3]==''and ir[5]=='':
                reactedA=intmed[ir[2]][2]*ir[0] #forward
                createdC=reactedA
                reactedC=intmed[ir[4]][2]*ir[1]#reverse
                createdA=reactedC
                #print(ir[8])
                
                #breakpoint()
                if ir[8].size>1 and ir[9].size>1:  #does this need to be remapped?
                    createdC = reactedA[ir[8]]
                    createdA = reactedC[ir[9]]
                if ir[2]=='citrate': #citrate -> isocitrate
                    if isocitr_channel==0: #two positions to add OH
                        createdC4 = 0.5*reactedA #Alex's isocitrate C4-OH
                        createdC2 = createdC4[ir[8]]
                        createdC=createdC2+createdC4
                        
                    elif isocitr_channel==2: #C2 for adding OH
                        createdC = reactedA[ir[8]]


                if ir[2]=='fumarate': #fumarate->malate
                    createdC2 = 0.5*reactedA
                    createdC3 = createdC2[ir[8]]
                    createdC=createdC2+createdC3

                intmed[ir[2]][3]+=createdA - reactedA
                intmed[ir[4]][3]+=createdC - reactedC
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)


        # Breakdown reaction  A->C+D
            elif ir[3]=='' and ir[5]!='':
                reactedA=intmed[ir[2]][2]*ir[0] #forward
                nC_C = intmed[ir[4]][0]

                createdC,createdD=breakdown(reactedA,nC_C)

                createdA,reactedC,reactedD=synthesis(intmed[ir[4]][2],intmed[ir[5]][2],ir[1])#reverse
                
                
                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    reactedM = reactedA[ir[8]]
                    createdC,createdD=breakdown(reactedM,nC_C)
                    createdA = createdA[ir[9]]

                intmed[ir[2]][3] +=createdA - reactedA
                intmed[ir[4]][3] +=createdC - reactedC
                intmed[ir[5]][3] +=createdD - reactedD
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)



        # Synthesis reaction  A+B->C
            elif ir[3]!='' and ir[5]=='':
                createdC,reactedA,reactedB=synthesis(intmed[ir[2]][2],intmed[ir[3]][2],ir[0])#forward

                reactedC = intmed[ir[4]][2]*ir[1]
                createdA,createdB=breakdown(reactedC,intmed[ir[2]][0])

                if ir[4]=='citrate': #acecoa+oxaace->citrate
                    createdC = createdC[ir[8]]
                    reactedC = reactedC[ir[9]]

                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    createdC = createdC[ir[8]]
                    reactedM = reactedC[ir[9]]
                    createdA,createdB=breakdown(reactedM,intmed[ir[2]][0])

                intmed[ir[2]][3] +=createdA -reactedA
                intmed[ir[3]][3] +=createdB - reactedB
                intmed[ir[4]][3] +=createdC - reactedC
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)
                
                if ir[2] == ir[3]:
                    flux_tracker[m][it] *=2

        # A + B --> C + D
            elif ir[2]!='' and ir[5]!='':
                createdM, reactedA,reactedB = synthesis(intmed[ir[2]][2],intmed[ir[3]][2],ir[0]) #make intermediate M
                rev_createdM,reactedC,reactedD = synthesis(intmed[ir[4]][2],intmed[ir[5]][2],ir[1]) #reverse reaction
                #breakpoint()
                if ir[8].size>1 and ir[9].size>1:    #does this need to be remapped?
                    createdM = createdM[ir[8]]
                    rev_createdM = rev_createdM[ir[9]]

                elif ir[5] == 'co2_als':
                    createdM = createdM[pyr_pyr_index]*0.5
                    rev_createdM = rev_createdM[pyr_pyr_index]

                nC_C=intmed[ir[4]][0]
                createdC,createdD=breakdown(createdM,nC_C)
                nC_A=intmed[ir[2]][0]
                createdA,createdB = breakdown(rev_createdM,nC_A)

                intmed[ir[2]][3] +=createdA - reactedA
                intmed[ir[3]][3] +=createdB - reactedB
                intmed[ir[4]][3] +=createdC - reactedC
                intmed[ir[5]][3] +=createdD - reactedD
                if it_check.is_integer():
                    flux_tracker[m][int(it_check)] = sum(reactedA)-sum(createdA)
                
                if ir[2] == ir[3] and it_check.is_integer():
                    flux_tracker[m][int(it_check)] *=2
                    
            m = m+1
    #changing concentration of intermediates
        i = 0
       
        for row in substrates:
            if sum(intmed[row][2]) > 0 and it_check.is_integer():
                concentration_tracker[i][int(it_check)]=sum(intmed[row][2])
                isotope_track = GetTotalRatio(intmed[row][0],intmed[row][2])
                isotope_tracker[i][int(it_check)] = GetDelta(isotope_track)
            intmed[row][2] += intmed[row][3]
            intmed[row][3] = np.zeros(2**intmed[row][0])
            if row in reservoirs:
                intmed[row][2][:] = intmed[row][4][:]
            
            i = i+1
        
    #const supply of the following reagents
            

        if PK_equilibration ==3:
            intmed['pyruvate'][2][:] = (intmed['pyruvate'][2][:]+intmed['PEP'][2][:])/2
            intmed['PEP'][2][:] = (intmed['pyruvate'][2][:]+intmed['PEP'][2][:])/2

        if it+2==t_steps:
            prev_intmed = deepcopy(intmed)
        
        flux_tracked = {}
            
    for i in range(len(reactions)):
        flux_tracked[reactions[i]] = flux_tracker[i][:]
    
    ## NETWORK DESCRIPTIONS:
    print('Number of reactions: ',len(EClist))     # Number of reactions
    print('Number of substrates: ',len(substrates)) # Number of substrates
    isotopologues = 0
    for row in substrates:
        isotopologues += 2**(intmed[row][0])
    print('Number of isotopologues: ', isotopologues)
    ## Print proportion of labelled substrates in 13C labelling simulations
    #for i in substrates:
       # M = GetProp(intmed[i][2])
       # print(i)
       # for j in np.arange(len(M)):
            #print('M',j,': ', M[j])
    return intmed, substrates, concentration_tracker, reactions, flux_tracker, metanetwork, EClist,flux_tracked, reservoirs, isotope_tracker, skip
def emprint(compound,intmed):
    value = siteSpecific(compound, intmed[compound][0], intmed[compound][2])
    return value
def prt3(str): # print out 3 effective digits
    print(str,'={0:.3f}'.format(eval(str)))
def prt10(str):
    print(str,'={0:.10f}'.format(eval(str)))
