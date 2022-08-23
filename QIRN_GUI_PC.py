from tkinter import *
from PIL import ImageTk,Image
from tkinter import filedialog
from Scripts import QIRN_20220218 as QIRNfile
from Scripts import QIRN_20211105_fluxinversion as inversionfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from Scripts import fluxcsveditor as fluxcsveditor

runcounter = 0
root = Tk()
root.title('Quantifying Isotopologue Reaction Networks (QIRN)')
#root.geometry('700x500')
root.configure(bg='white')
def runinversion():
    global frame
    global my_img
    global intmed
    global substrates
    global time_steps
    global concentration_tracker
    global filename
    global time
    global TimeEntry
    global reactions
    global flux_button
    global concentration_button
    global isotope_button
    global again_button
    global runcounter
    
    if runcounter > 0:
        flux_button.destroy()
        isotope_button.destroy()
        concentration_button.destroy()
        again_button.destroy()
    
    filename = filedialog.askopenfilename(initialdir = "/Documents/QIRN/Photosynthesis", title = "Select Network File")
    TimeEntry = Entry(root)
    TimeEntry.grid(row = 2, column = 1,pady=5, columnspan = 3 )
    TimeEntry.insert(0,"Enter Model Time (sec):")
    
    
    runButton = Button(root,text="Invert Fluxes", command=invertfluxes)
    runButton.grid(row=3,column=0, pady=5)


def runQIRN():
    global frame
    global my_img
    global intmed
    global substrates
    global time_steps
    global concentration_tracker
    global filename
    global time
    global TimeEntry
    global runButton
    global reactions
    global flux_button
    global concentration_button
    global isotope_button
    global isotope_time_button
    global again_button
    global runcounter
    if runcounter > 0:
        flux_button.destroy()
        isotope_button.destroy()
        concentration_button.destroy()
        again_button.destroy()
        isotope_time_button.destroy()
    
    filename = filedialog.askopenfilename(initialdir = "/Documents/QIRN/Photosynthesis", title = "Select a Network File")
    TimeEntry = Entry(root)
    TimeEntry.grid(row = 2, column = 1,pady=5, columnspan=2)
    TimeEntry.insert(0,"Enter Model Time (sec):")
    
    
    runButton = Button(root,text="Run Model", command=run)
    runButton.grid(row=3,column=1, pady=5, columnspan =2)


def invertfluxes():
    global filename
    global frame
    global my_img
    global intmed
    global substrates
    global time_steps
    global concentration_tracker
    global time
    global TimeEntry
    global metanetwork
    global reactions
    global flux_tracker
    global dt
    global flux_button
    global concentration_button
    global isotope_button
    global again_button
    global fluxbutton_photo
    global concbutton_photo
    global isotopebutton_photo
    global fluxtracking
    global EClist
    global reservoirs
    global K
    global inversionfile
    global kr
    global kf
    time = int(TimeEntry.get())
    dt = 0.1
    kf, kr = inversionfile.fluxinvert("IntermediatesDatabase.csv",filename,"ReactionDatabase.csv",time,dt)
    print(kf)
    print(kr)
    fluxcsveditor.fluxconvert(kf,kr,filename)
    
def run():
    global filename
    global frame
    global my_img
    global intmed
    global substrates
    global time_steps
    global concentration_tracker
    global time
    global TimeEntry
    global metanetwork
    global reactions
    global flux_tracker
    global dt
    global flux_button
    global concentration_button
    global isotope_button
    global isotope_time_button
    global again_button
    global fluxbutton_photo
    global concbutton_photo
    global isotopebutton_photo
    global runButton
    global isotopeTIMEbutton_photo
    global fluxtracking
    global EClist
    global reservoirs
    global runcounter
    global K
    global kr
    global kf
    global isotope_tracker
    global skip

    time = int(TimeEntry.get())
    dt = 0.1
    intmed,substrates,concentration_tracker, reactions, flux_tracker, metanetwork, EClist,fluxtracking, reservoirs, isotope_tracker,skip = QIRNfile.QIRN("IntermediatesDatabase.csv",filename,"ReactionDatabase.csv",time,dt)
    time_steps = int(time/dt)
    fluxbutton_photo = Image.open('GUI_Illustrations/reaction.png')
    fluxbutton_photo = fluxbutton_photo.resize((200,35), Image.ANTIALIAS)
    fluxbutton_photo = ImageTk.PhotoImage(fluxbutton_photo)
    concbutton_photo = Image.open('GUI_Illustrations/intermediate.png')
    concbutton_photo = concbutton_photo.resize((200,35), Image.ANTIALIAS)
    concbutton_photo = ImageTk.PhotoImage(concbutton_photo)
    
    isotopebutton_photo = Image.open('GUI_Illustrations/isotope.png')
    isotopebutton_photo = isotopebutton_photo.resize((200,35), Image.ANTIALIAS)
    isotopebutton_photo = ImageTk.PhotoImage(isotopebutton_photo)
    
    isotopeTIMEbutton_photo = Image.open('GUI_Illustrations/CSIAwTIME.png')
    isotopeTIMEbutton_photo = isotopeTIMEbutton_photo.resize((200,35), Image.ANTIALIAS)
    isotopeTIMEbutton_photo = ImageTk.PhotoImage(isotopeTIMEbutton_photo)
    
    
    isotope_button = Button(root, image = isotopebutton_photo, command = lambda: open('isotopes'),borderwidth=0,bg='white',highlightthickness=0,bd=0)
    isotope_button.grid(row=3,column=2)
    
    isotope_time_button = Button(root, image = isotopeTIMEbutton_photo, command = lambda: open('isotope_time'),borderwidth=0,bg='white',highlightthickness=0,bd=0)
    isotope_time_button.grid(row=3,column=3)
   
    flux_button = Button(root, image = fluxbutton_photo, command = lambda: open('fluxes'),borderwidth=0,bg='white',highlightthickness=0,bd=0)
    flux_button.grid(row=3,column=1)
    
    concentration_button = Button(root, image = concbutton_photo, command = lambda: open('concentrations'),borderwidth=0,bg='white',highlightthickness=0,bd=0)
    concentration_button.grid(row=3,column=0)
    
    again_button = Button(root, text = "Run Again", command = runagain)
    again_button.grid(row=4,column=1, pady=5, columnspan =2)

    runButton.destroy()

def runagain():
    global filename
    global frame
    global my_img
    global intmed
    global substrates
    global time_steps
    global concentration_tracker
    global time
    global TimeEntry
    global metanetwork
    global reactions
    global runcounter
    global flux_tracker
    global dt
    global flux_button
    global concentration_button
    global isotope_button
    global again_button
    global isotope_time_button
    global fluxbutton_photo
    global fluxbutton_photo
    global concbutton_photo
    global isotopebutton_photo
    global isotopeTIMEbutton_photo
    global isotope_time_button
    global runButton
    global fluxtracking
    global EClist
    global reservoirs
    global kf
    global kr
    global isotope_tracker
    
    flux_button.destroy()
    runButton.destroy()
    isotope_button.destroy()
    concentration_button.destroy()
    isotope_time_button.destroy()
    time = int(TimeEntry.get())
    dt = 0.1
    intmed,substrates,concentration_tracker, reactions, flux_tracker, metanetwork, EClist,fluxtracking, reservoirs, isotope_tracker,skip = QIRNfile.QIRN("IntermediatesDatabase.csv",filename,"ReactionDatabase.csv",time,dt)
    time_steps = int(time/dt)

    isotope_button = Button(root, image = isotopebutton_photo, command = lambda: open('isotopes'))
    isotope_button.grid(row=3,column=2)
    flux_button = Button(root, image = fluxbutton_photo, command = lambda: open('fluxes'))
    flux_button.grid(row=3,column=1)
    
    isotope_time_button = Button(root, image = isotopeTIMEbutton_photo, command = lambda: open('isotope_time'))
    isotope_time_button.grid(row=3,column=3)
   
    concentration_button = Button(root, image = concbutton_photo, command = lambda: open('concentrations'))
    concentration_button.grid(row=3,column=0)
    runcounter = runcounter+1

row_num =1

def clear():
    global Isotopes
    Isotopes.destroy()
    Isotopes = LabelFrame(top,padx=20)
    Isotopes.grid(row = 0, column=3)
def clicked(molecule):
    global intmed
    global top
    global Intermediates
    global logo2
    global Isotopes
    global row_num
    global MODES
    global reactions

    nC = intmed[molecule][0]
    conc = intmed[molecule][2] 
    
    PSIALabel = Label(Isotopes,text = molecule,bg='white')
    PSIALabel.grid(column = 3, row =row_num)
    row_num= row_num+1
    for iC in range(nC):
        cumu_ratio=QIRNfile.GetSingleRatio(conc,iC)
        cumu_dC=QIRNfile.GetDelta(cumu_ratio)
        #print("C{}: {} ".format(iC+1, cumu_dC))
        PSIALabel = Label(Isotopes,text=("C{}: {} ".format(iC+1, cumu_dC)),bg='white')
        PSIALabel.grid(column = 3, row =row_num)
        row_num=row_num+1
fluxcounter = 0

def run_fluxes_plot():
    global MODES
    global substrates
    global vars_flux
    global time_steps
    global concentration_tracker
    global Flux
    global fig
    global metanetwork
    global reactions
    global fluxframe
    global flux_tracker
    global dt
    global fluxtracking
    global fluxcounter
    global skip
    
    flux_plot_indices =[]
    flux_plot_reactions =[]
    flux_plot_columns = []
   
    for i in range(len(reactions)):
        if vars_flux[i].get() == 1:
            flux_plot_indices.append(i)
    for i in flux_plot_indices:
        flux_plot_reactions.append(reactions[i])
    for i in range(len(reactions)):
        if reactions[i] in flux_plot_reactions:
            flux_plot_columns.append(i)
    #breakpoint()
    time_window = np.arange(0,time_steps+1,skip)
    flux_tracked = flux_tracker.transpose()
    hfont = {'fontname':'Helvetica','fontweight':'bold'}
    
    fig = plt.figure()
    fig = Figure(figsize=(6, 6), dpi=100)
    #fig.add_subplot().plot(concentration_plotted)
    ax = fig.add_subplot(111)
    rownum = 1
    for i in flux_plot_columns:
        #breakpoint()
        ax.plot(time_window*dt,flux_tracked[:,i]*(1/dt),label = str(reactions[i]))
        endflux = "{:.4f}".format(fluxtracking[reactions[i]][-1]*100)
        text = str(reactions[i])+" = " +str(endflux)
        fluxlabel = Label(fluxframe, text = text)
        fluxlabel.grid(row = rownum, column =1)
        rownum = rownum+1
    ax.legend()
    ax.set_xlabel('Time (sec)',fontname = 'Helvetica',fontweight = 'bold')
    ax.set_ylabel('Flux (mol/sec)', **hfont)
    ax.set_title('Reaction Fluxes',**hfont)
    canvas = FigureCanvasTkAgg(fig, master = fluxframe)
    canvas.draw()
    canvas.get_tk_widget().grid(column = 1, row =0,padx = 5)
    
def run_fluxes(var1,var2):
    global MODES
    global substrates
    global vars
    global time_steps
    global concentration_tracker
    global isotope_tracker
    global Concentration
    global fig1
    global plotframe
    global reactions
    global top
    global isotope_tracker
    global top1
    global skip
    concentration_plot_indices =[]
    concentration_plot_substrates =[]
    concentration_plot_columns = []
    plotframe.destroy()
    plotframe = LabelFrame(top1, padx=150,bg='white',bd=0)
    plotframe.grid(row = 0, column = 2, rowspan = 10, padx = 60, pady=5)
    #plotframe = LabelFrame(top, padx=20)
    #plotframe.grid(row = 0, column = 2, rowspan = 10, padx = 20, pady=20)
    for i in range(len(MODES)):
        if vars[i].get() == 1:
            concentration_plot_indices.append(i)
    for i in concentration_plot_indices:
        concentration_plot_substrates.append(MODES[i])
    for i in range(len(substrates)):
        if substrates[i] in concentration_plot_substrates:
            concentration_plot_columns.append(i)
    
    time_window = np.arange(0,time_steps+1,skip)
    concentration_tracked = var1.transpose()
    #concentration_data = pd.DataFrame(concentration_tracked[time_window,:]*(1/dt),index = time_window*dt,columns = substrates)
    
    #concentration_plotted = concentration_data.iloc[:,concentration_plot_columns]
    #MProportions_biosynth = MProportions_data.iloc[:,58:-1]
    hfont = {'fontname':'Helvetica','fontweight':'bold'}
    
    fig1 = plt.figure()
    fig1 = Figure(figsize=(6, 6), dpi=100)
    #fig.add_subplot().plot(concentration_plotted)
    ax1 = fig1.add_subplot(111)

    rownum = 1
    for i in concentration_plot_columns:
        ax1.plot(time_window*dt,concentration_tracked[:,i],label = substrates[i])
        endflux = "{:.1f}".format(concentration_tracked[-1][i])
        text = str(substrates[i])+ " = " +str(endflux)
        fluxlabel = Label(plotframe, text = text)
        fluxlabel.grid(row = rownum, column =1)
        rownum = rownum+1
        
    ax1.legend()
    ax1.set_xlabel('Time Steps',fontname = 'Helvetica',fontweight = 'bold')
    if var2 == 'conc':
        ax1.set_ylabel('Concentration (mol/volume)', **hfont)
        ax1.set_title('Concentration of Substrates',**hfont)
    elif var2 == 'isotope':
        ax1.set_ylabel('Isotope Composition (‰)', **hfont)
        ax1.set_title('Compund-Specific Isotope Composition',**hfont)

    canvas1 = FigureCanvasTkAgg(fig1, master = plotframe)
    canvas1.draw()
    canvas1.get_tk_widget().grid(column = 1, row =0,padx = 5)
    
    
def check(name):
    print(name)
def open(window):
    global intmed
    global top
    global substrates
    global frame
    global my_img1
    global Concentration
    global Intermediates
    global logo2
    global Isotopes
    global MODES
    global vars
    global vars_flux
    global Flux
    global plotframe
    global reactions
    global fluxframe
    global isotope_tracker
    global concentration_tracker
    global top1

    if window == 'isotopes':
        top = Toplevel(bg = 'white')
        top.title('Isotope Distribution of Intermediates')
        Intermediates = LabelFrame(top, padx=20,bg='white',bd=0)
        Intermediates.grid(row = 0, column = 0, rowspan = 10)
        Isotopes = LabelFrame(top,padx=20,bg='white')
        Isotopes.grid(row = 0, column=3)
        
        MODES =[]

        for i in substrates:
            if intmed[i][2][0] > 0:
                MODES.append(i)
        molecule = StringVar()
        j=0
        r = 0
        for text in MODES:
            Radiobutton(Intermediates,text=text, variable = molecule, value = text,bg='white').grid(column=int(r), row=j)
            j=j+1
          
            if j ==15:
                r += 1
                j = 0
        PSIA_button = Button(Intermediates,text="Calculate Isotope Distribution", bg = 'white', command=lambda: clicked(molecule.get()))
        PSIA_button.grid(column=int(r/2),row=16,pady=10)
        clear_button = Button(Intermediates,text="Clear", command=clear, bg = 'white')
        clear_button.grid(column=int(r/2),row=17,pady=10)
        
    if window == 'concentrations':
        top1 = Toplevel(bg='white')
        top1.title('Concentration of Intermediates')
        Concentration = LabelFrame(top1, padx=20,bg='white',bd=0)
        Concentration.grid(row = 0, column = 0, rowspan = 10, padx = 20, pady=20)
        plotframe = LabelFrame(top1, padx=60,bg='white',bd=0)
        plotframe.grid(row = 0, column = 2, rowspan = 10, padx = 60, pady=20)
        MODES =[]

        for i in substrates:
            if intmed[i][2][0] > 0:
                MODES.append(i)
        
        j=0
        r = 0
        vars = []
        for text1 in MODES:
            var = IntVar()
            vars.append(var)
            Checkbutton(Concentration,text=text1, variable = var,bg='white').grid(column=r, row=j)
            j=j+1
            if j ==15:
                r += 1
                j = 0
        concentration_button = Button(Concentration,text="Display Concentrations with Time", command = lambda: run_fluxes(concentration_tracker,'conc'))
        concentration_button.grid(column=int(r/2),row=17, pady=10)
    
    if window == 'isotope_time':
        top1 = Toplevel(bg='white')
        top1.title('Molecular Average Isotope Composition of Intermediates')
        Concentration = LabelFrame(top1, padx=20,bg='white',bd=0)
        Concentration.grid(row = 0, column = 0, rowspan = 10, padx = 20, pady=20)
        plotframe = LabelFrame(top1, padx=60,bg='white',bd=0)
        plotframe.grid(row = 0, column = 2, rowspan = 10, padx = 60, pady=20)

        MODES =[]

        for i in substrates:
            if intmed[i][2][0] > 0:
                MODES.append(i)
        
        j=0
        r = 0
        vars = []
        for text1 in MODES:
            var = IntVar()
            vars.append(var)
            Checkbutton(Concentration,text=text1, variable = var,bg='white').grid(column=int(r), row=j)
            j=j+1
            if j ==15:
                r += 1
                j = 0
        concentration_button = Button(Concentration,text="Display Isotope Composition with Time", command = lambda: run_fluxes(isotope_tracker,'isotope'))
        concentration_button.grid(column=int(r/2),row=17,pady=10)
    
    if window == 'fluxes':
        top2 = Toplevel(bg='white')
        top2.title('Flux of Reactions')
        Flux = LabelFrame(top2, padx=20,bg='white',bd=0)
        Flux.grid(row = 0, column = 0, rowspan = 10, padx = 5, pady=20)
        fluxframe = LabelFrame(top2, padx=150,bg='white',bd=0)
        fluxframe.grid(row = 0, column = 2, rowspan = 10, padx = 10, pady=20)
        logoframe = LabelFrame(top2, padx = 60, bd = 0, relief = SUNKEN,bg='white')
        logoframe.grid(row = 11, column = 0)
        #MODES =[]
        j=0
        r = 0
        vars_flux = []
        for text in reactions:
            var1 = IntVar()
            vars_flux.append(var1)
            Checkbutton(Flux,text=text, variable = var1,bg='white').grid(column=r, row=j)
            j=j+1
            if j ==15:
                r += 1
                j = 0
        flux_button = Button(Flux,text="Display Reaction Fluxes with Time", command=run_fluxes_plot)
        flux_button.grid(column=int(r/2),row=17,pady=10)
        
        logo1 = Image.open('GUI_Illustrations/bluelogo.png')
        logo1 = logo1.resize((300,150), Image.ANTIALIAS)
        my_img1 = ImageTk.PhotoImage(logo1)
        logo1 = Label(logoframe, image = my_img1,bg='white').grid(row =0,column=0)
        
        
#frame = LabelFrame(root, padx=10,pady=10)
#frame.pack(padx=5,pady=5)


logo = Image.open('GUI_Illustrations/QIRN.png')
logo = logo.resize((700,575), Image.ANTIALIAS)
my_img = ImageTk.PhotoImage(logo)
logo = Label(root, image = my_img,bg='white').grid(row=0,column=0,columnspan=4)
button_photo = Image.open('GUI_Illustrations/submitbutton3.png')
button_photo = button_photo.resize((220,35), Image.ANTIALIAS)
button_photo = ImageTk.PhotoImage(button_photo)
ConvertButton = Image.open('GUI_Illustrations/ConvertButton.png')
ConvertButton = ConvertButton.resize((240,34), Image.ANTIALIAS)
ConvertButton = ImageTk.PhotoImage(ConvertButton)

submit_button = Button(root, image = button_photo, command = runQIRN,borderwidth=0,bg='white',highlightthickness=0,bd=0)
submit_button.grid(row = 1, column=2, padx = 0, pady=10,columnspan=1)

inversion_button = Button(root, image = ConvertButton, command = runinversion,borderwidth=0,bg='white',highlightthickness=0,bd=0)
inversion_button.grid(row = 1, column=1, padx = 0, pady=10,columnspan=1)
root.mainloop()





