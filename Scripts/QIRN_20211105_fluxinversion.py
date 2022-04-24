#Elliott Mueller, GPS Caltech, 2022
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
        print(singlerates)
        return np.repeat(singlerates,2**nC)
    else:
        print(singlerates)
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
    
def fluxinvert(intermediates,networkbuilder,reactiondatabase, time, dt):

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

    with open(intermediates, newline='') as imfile:
            imreader = csv.reader(imfile)
            next(imreader)  #get rid of the first row
            for row in imreader:
                substrates.append(row[0])


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
  

  
    #INVERT MASS ACTION RATE CONSTANTS FROM ABSOLUTE FLUXES
    #OCTOBER 21 2021


    counter = 0
    EC_list = EClist

    requestedflux = np.zeros(len(EC_list))
    rev_requestedflux = np.zeros(len(EC_list))
    reversibility = np.zeros(len(EC_list))
    requestedreversibility = np.zeros(len(EC_list))

    #Kf = np.ones(len(EC_list)*2)
    #Putting reverse reaction rates below forward ones on the same vector... You can split vector in half always to separate.
    
    time_steps = int(time)
    ts = np.linspace(0, time, time_steps)

    #reservoirs = [0]
    substrates =[]
    reservoir_index = []
    kindex = []
    
    for row in EC_list:
        ir = metanetwork[row]
        if ir[2] not in substrates and ir[2] != '':
            substrates.append(ir[2])
        if ir[3] not in substrates and ir[3] != '':
            substrates.append(ir[3])
        if ir[4] not in substrates and ir[4] != '':
            substrates.append(ir[4])
        if ir[5] not in substrates and ir[5] != '':
            substrates.append(ir[5])
    
    conc0 = np.zeros(len(substrates))
    conc = np.zeros(len(substrates))
    
    rxnflux = np.zeros([len(substrates),len(EC_list)])
    rev_rxnflux = np.zeros([len(substrates),len(EC_list)])
    rev_reactantA = np.zeros([len(substrates),len(EC_list)])
    rev_reactantB = np.ones([len(substrates),len(EC_list)])
    reactantA = np.zeros([len(substrates),len(EC_list)])
    reactantB = np.ones([len(substrates),len(EC_list)])
    rev_rxnindex = np.zeros([len(substrates),len(EC_list)])
    rxnindex = np.zeros([len(substrates),len(EC_list)])
    rev_secondaryrxnindex = np.zeros([len(substrates),len(EC_list)])
    secondaryrxnindex = np.zeros([len(substrates),len(EC_list)])
    
    
    for row in substrates:
        conc0[substrates.index(row)] = intmed[row][1]
        conc[substrates.index(row)] = intmed[row][1]

    for i in reservoirs:
        if i in substrates:
            reservoir_index = np.append(reservoir_index,substrates.index(i))
            print(reservoir_index)

    for row in EC_list:
        ir = metanetwork[row]
        requestedflux[counter] = float(ir[0])
        rev_requestedflux[counter] = float(ir[1])
        A = substrates.index(ir[2])
        C = substrates.index(ir[4])
        
        rxnflux[A,counter] = -1
        rxnflux[C,counter] = 1
        rev_rxnflux[A,counter] = -1
        rev_rxnflux[C,counter] = 1

        if ir[3] != '':
            B = substrates.index(ir[3])
            rxnflux[B,counter] = -1
            rev_rxnflux[B,counter] = -1
                
            if A == B:
                rxnflux[A,counter] = -2
                rev_rxnflux[A,counter] = -2
       
        if ir[5] != '':
            D = substrates.index(ir[5])
            rxnflux[D,counter] = 1
            rev_rxnflux[D,counter] = 1
            
            if C == D:
                rxnflux[C,counter] = 2
                rev_rxnflux[C,counter] = 2
        
        counter += 1


    print(substrates)
    
    for i in range(len(EClist)):
        if requestedflux[i] != 0 and rev_requestedflux[i] != 0:
            requestedreversibility[i] = rev_requestedflux[i]/requestedflux[i]
        if requestedflux[i] == 0:
            rxnflux[:,i] = 0
        else:
            kindex = np.append(kindex,i)
    
    for i in range(len(EClist)):
        if rev_requestedflux[i] == 0:
            rev_rxnflux[:,i] = 0
        else:
            kindex = np.append(kindex,(i+int(len(EClist))))
    kstart = []
    
    count = 0
    for row in EC_list:
        ir = metanetwork[row]
        if ir[0] != '0':
            kstart = np.append(kstart, float(ir[0]))
        count += 1
        
    count = 0
    for row in EC_list:
        ir = metanetwork[row]
        if ir[1] != '0':
            kstart = np.append(kstart, float(ir[1]))
        count += 1
   
    Kf = kstart
    Kconst = (Kf)
    tracker = np.zeros(len(substrates))
    print(rxnflux)
    print('')
    print(rev_rxnflux)
    print('')
    print(kindex)
    print(kstart)
    
    for j in range(len(substrates)):
        for i in range(len(EClist)):
            if rxnflux[j,i] == 1:
                for k in range(len(rxnflux[:,i])):
                    if rxnflux[k,i] == -1 and j != k:
                        rxnindex[j,i] = k+1
                for m in range(len(rxnflux[:,i])):
                    if rxnflux[m,i] == -1 and rxnindex[j,i] != m+1:
                        secondaryrxnindex[j,i] = rxnindex[j,i]
                        rxnindex[j,i] = m+1
                    if rxnflux[m,i] == -2 and rxnindex[j,i] != m+1:
                        secondaryrxnindex[j,i] = m+1
                        rxnindex[j,i] = m+1
                        
            elif rxnflux[j,i] == -1:
                rxnindex[j,i] = j+1
                for k in range(len(rxnflux[:,i])):
                    if rxnflux[k,i] == -1 and j != k:
                        secondaryrxnindex[j,i] = k+1
            
            
            elif rxnflux[j,i] == 2:
                for k in range(len(rxnflux[:,i])):
                    if rxnflux[k,i] == -1:
                        rxnindex[j,i] = k+1
                    elif rxnflux[k,i] == -2:
                        rxnindex[j,i] = k+1
                        secondaryrxnindex[j,i] = k+1
                for m in range(len(rxnflux[:,i])):
                    if rxnflux[m,i] == -1 and rxnindex[j,i] != m+1:
                        secondaryrxnindex[j,i] = rxnindex[j,i]
                        rxnindex[j,i] = m+1
                        
                        
            elif rxnflux[j,i] == -2:
                rxnindex[j,i] = j+1
                secondaryrxnindex[j,i] = j+1
            
            if rev_rxnflux[j,i] == 1:
                rev_rxnindex[j,i] = j+1
                for k in range(len(rxnflux[:,i])):
                    if rev_rxnflux[k,i] == 1 and j != k:
                        rev_secondaryrxnindex[j,i] = k+1
                        
            elif rev_rxnflux[j,i] == -1:
                for k in range(len(rxnflux[:,i])):
                    if rev_rxnflux[k,i] == 1 and j != k:
                        rev_rxnindex[j,i] = k+1
                for m in range(len(rxnflux[:,i])):
                    if rev_rxnflux[m,i] == 1 and rev_rxnindex[j,i] != m+1:
                        rev_secondaryrxnindex[j,i] = rev_rxnindex[j,i]
                        rev_rxnindex[j,i] = m+1
                    if rev_rxnflux[m,i] == -2 and rev_rxnindex[j,i] != m+1:
                        rev_secondaryrxnindex[j,i] = rev_rxnindex[j,i]
                        rxnindex[j,i] = m+1
            
            elif rev_rxnflux[j,i] == -2:
                for k in range(len(rxnflux[:,i])):
                    if rev_rxnflux[k,i] == 1:
                        rev_rxnindex[j,i] = k+1
                    elif rev_rxnflux[k,i] == 2:
                        rev_rxnindex[j,i] = k+1
                        rev_secondaryrxnindex[j,i] = k+1
                for m in range(len(rxnflux[:,i])):
                    if rev_rxnflux[m,i] == 1 and rev_rxnindex[j,i] != m+1:
                        rev_secondaryrxnindex[j,i] = rev_rxnindex[j,i]
                        rev_rxnindex[j,i] = m+1
                        
            elif rev_rxnflux[j,i] == 2:
                rev_rxnindex[j,i] = j+1
                rev_secondaryrxnindex[j,i] = j+1
    
    #rxnflux = rxnflux.T
    #rev_rxnflux = rev_rxnflux.T
    #rev_reactantA = rev_reactantA.T
    #rev_reactantB = rev_reactantB.T
    #reactantA = reactantA.T
    #reactantB = reactantB.T
    #rev_rxnindex = rev_rxnindex.T
    #rxnindex = rxnindex.T
    #rev_secondaryrxnindex = rev_secondaryrxnindex.T
    #secondaryrxnindex = secondaryrxnindex.T
    
    def dK_dt(conc, t,Kf):
        Kvector = np.zeros(len(EClist)*2)
        count = 0
        for i in kindex:
            Kvector[int(i)] = Kf[count]
            count += 1
            
        Kf = Kvector
        Kr = Kf[int(len(Kf)/2)::]
        Kf = Kf[0:int(len(Kf)/2)]
    
        for i in range(len(rxnindex[:,0])):
            for j in range(len(rxnindex[0,:])):
                if rxnindex[i,j] != 0:
                    reactantA[i,j] = conc[int(rxnindex[i,j])-1]
                if secondaryrxnindex[i,j] != 0:
                    reactantB[i,j] = conc[int(secondaryrxnindex[i,j])-1]
                if rev_rxnindex[i,j] != 0:
                    rev_reactantA[i,j] = conc[int(rev_rxnindex[i,j])-1]
                if rev_secondaryrxnindex[i,j] != 0:
                    rev_reactantB[i,j] = conc[int(rev_secondaryrxnindex[i,j])-1]
        tracker = np.matmul((rxnflux * reactantA * reactantB),Kf) + np.matmul((-1*rev_rxnflux * rev_reactantA * rev_reactantB),Kr)
        for n in reservoir_index:
            tracker[int(n)] = 0
        
        return tracker

    def f(Kconst):
        K = odeint(dK_dt, conc0, ts, args = (Kconst,))
        
        
        flux = np.zeros(len(EClist))
        rev_flux = np.zeros(len(EClist))
        Kvector = np.zeros(len(EClist)*2)
        past_flux = np.zeros(len(EClist))
        past_rev_flux = np.zeros(len(EClist))
        past_rev_flux = np.zeros(len(EClist))
        
        Kvector = np.zeros(len(EClist)*2)
        count = 0
        for i in kindex:
            Kvector[int(i)] = Kconst[count]
            count += 1
            
        Kconst = Kvector
        
        Kconst_r = Kconst[int(len(Kconst)/2)::]
        Kconst = Kconst[0:int(len(Kconst)/2)]


        for j in range(len(substrates)):
            for i in range(len(EClist)):
                intermediate = 0
                rev_intermediate = 0
                past_intermediate = 0
                revpast_intermediate = 0
               # if flux[i] == 0:
                    
                if rxnflux[j,i] == -1:
                    intermediate = Kconst[i]
                    intermediate *= K[-1,j]
                    past_intermediate = abs(Kconst[i]*K[-100,j])

                    for k in range(len(rxnflux[:,i])):
                        if rxnflux[k,i] == -1 and j != k:
                            intermediate *= K[-1,k]
                            past_intermediate *= abs(K[-100,k])
                    flux[i] = abs(intermediate)
                    past_flux[i] = abs(past_intermediate)
                    
                elif rxnflux[j,i] == -2:
                    intermediate = Kconst[i]
                    intermediate *= (K[-1,j]**2)
                    flux[i] = abs(intermediate)
                    past_flux[i] = abs((K[-100,j])**2) * Kconst[i]
            
           # if rev_flux[i] == 0:
                if rev_rxnflux[j,i] == 1:
                    rev_intermediate = Kconst_r[i]
                    rev_intermediate *= K[-1,j]
                    past_revintermediate = abs(K[-100,j])*Kconst_r[i]
                    
                    for k in range(len(rxnflux[:,i])):
                        if rev_rxnflux[k,i] == 1 and j != k:
                            rev_intermediate *= K[-1,k]
                            past_revintermediate *= abs(K[-100,k])

                    rev_flux[i] = abs(rev_intermediate)
                    past_rev_flux[i] = abs(past_revintermediate)
                    
                elif rev_rxnflux[j,i] == 2:
                    rev_intermediate = Kconst_r[i]
                    rev_intermediate *= K[-1,j]**2
                    rev_flux[i] = abs(rev_intermediate)
                    past_rev_flux[i] = (K[-100,j]**2)*Kconst_r[i]
                    
                    
                if rev_flux[i] > 1e-5 and flux[i] > 1e-5:
                    reversibility[i] = rev_flux[i]/flux[i]
                    if reversibility[i] > 1e3:
                        reversibility[i] = 0


        print('')

        for i in range(len(flux)):
            if float(flux[i]) == 0:
                flux[i] = 1
            if float(rev_flux[i]) == 0:
                rev_flux[i] = 1
            if float(requestedflux[i]) == 0:
                requestedflux[i] = 1
            if float(rev_requestedflux[i]) == 0:
                rev_requestedflux[i] = 1
            if float(past_rev_flux[i]) == 0:
                past_rev_flux[i] = 1
            if float(past_flux[i]) == 0:
                past_flux[i] = 1
            if float(reversibility[i]) == 0:
                reversibility[i] = 1
                print(i,' HERE')
            if float(requestedreversibility[i]) == 0:
                requestedreversibility[i] = 1
    


        print('')
        print(flux/requestedflux)
        print(rev_flux/rev_requestedflux)
        print(flux/past_flux)
        print(rev_flux/past_rev_flux)
        print(reversibility/requestedreversibility)
        print(reversibility)
        print(requestedreversibility)
       # breakpoint()
        residual = sum(((flux)/(requestedflux)-1)**2)**0.5 + sum(((rev_flux)/(rev_requestedflux)-1)**2)**0.5 +  + 0.5*sum(((flux)/(past_flux)-1)**2)**0.5 + 0.5*sum(((rev_flux)/(past_rev_flux)-1)**2)**0.5 + sum(((reversibility)/(requestedreversibility)-1)**2)**0.5
        return residual
    
    
    bnd = (0,1)
    bnds = ((bnd,)*len(kindex))
    minimum = optimize.minimize(f,Kconst,method='SLSQP', bounds = bnds, tol = 1e-6)

    

    if minimum.fun > 0.1:
        print("No steady state found with those fluxes")
        print("Steady state found at fluxes  ", minimum.fun, "off from the specified values")
    else:
        print("Steady state found with those fluxes!")
        print("QIRN matched your specified fluxes by: ", minimum.fun)
    
    optvalindex = minimum.x
    optval = np.zeros(len(EClist)*2)
   # optval_f = np.array([0, 0, 0.1, 0, 0.0046, 0, 0, 0, 0, 0, 0, 0.05, 0, 0.05, 0, 0, 0.0092, 0, 0.0046, 0.0056, 0.0138, 0.0056, 0.025])
    
  #  optval_r = np.array([0.025, 0.0556, 0, 0.275, 0, 0.05, 0.05, 0.025, 0.0556, 0.275, 0.025, 0, 0.0046,0, 0.1, 0.1194, 0, 0.05, 0, 0, 0, 0, 0])
   # optvalindex = np.concatenate([optval_f,optval_r])
    K = odeint(dK_dt, conc0, ts, args = (optvalindex,))
    
    count = 0
    for i in kindex:
        optval[int(i)] = optvalindex[count]
        count += 1
        
    optval_r = optval[int(len(optval)/2)::]
    optval_f = optval[0:int(len(optval)/2)]
    for i in range(len(EClist)):
        if rev_requestedflux[i] == 0:
            optval_r[i] = 0
        if requestedflux[i] == 0:
            optval_f[i] = 0
    

    print(optval_f)
    print(optval_r)
    
    return optval_f,optval_r


