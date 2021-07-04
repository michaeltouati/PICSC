# Particle-In-Cell (PIC) Simulation Cost Estimator 
# written by Dr Michael J TOUATI - 06/28/2019 
# mtouati@clpu.es

import tkinter
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
import os

def estimate(*args):
	try:
		Ndim_value         = str(Ndim.get())
		Nspecies_value     = float(Nspecies.get())
		Npart_value        = float(Npart.get())
		Nt_value           = float(Nt.get())
		Interp_order_value = str(Interp_order.get())
		Nx_value           = float(Nx.get())
		Ny_value           = float(Ny.get())
		Nz_value           = float(Nz.get())
		# Evaluate the needed RAM memory and the best memory unit to express it
		# assuming double-floating-point format (1 number needs 64 bites or 8 B)
		doubles = 8. 
		if Ndim_value == '1D':
			res_part   = 4.*Npart_value*doubles                         # 4 for x, px, py and pz
			res_fields = (1.+((Nspecies_value+1.)*2.))*Nx_value*doubles # 3 for Ex, jx and rho
		elif Ndim_value == '2D':
			res_part   = 5.*Npart_value*doubles                                  # 5 for x, y, px, py and pz
			res_fields = (3.+((Nspecies_value+1.)*3.))*Nx_value*Ny_value*doubles # 6 for Ex, Ey, Bz, jx, jy and rho
		else:
			res_part   = 6.*Npart_value*doubles                                           #  6 for x, y, z, px, py and pz
			res_fields = (6.+((Nspecies_value+1.)*4.))*Nx_value*Ny_value*Nz_value*doubles # 10 for Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz and rho
		if res_part<1.e3:
			units_part = "B"
		elif res_part>=1.e3 and res_part<1.e6:
			units_part = "KB"
			res_part = res_part*1.e-3
		elif res_part>=1.e6 and res_part<1.e9:
			units_part = "MB"
			res_part = res_part*1.e-6
		elif res_part>=1.e9 and res_part<1.e12:
			units_part = "GB"
			res_part = res_part*1.e-9
		else:
			units_part = "TB"
			res_part = res_part*1.e-12
		if res_fields<1.e3:
			units_fields = "B"
		elif res_fields>=1.e3 and res_fields<1.e6:
			units_fields = "KB"
			res_fields = res_fields*1.e-3
		elif res_fields>=1.e6 and res_fields<1.e9:
			units_fields = "MB"
			res_fields = res_fields*1.e-6
		elif res_fields>=1.e9 and res_fields<1.e12:
			units_fields = "GB"
			res_fields = res_fields*1.e-9
		else:
			units_fields = "TB"
			res_fields = res_fields*1.e-12
		# experimental correction factors using one particular code
		# keeping only 3 decimals
		res_part   = int(1.e3*1.2537*res_part)/1.e3
		res_fields = int(1.e3*1.9278*res_fields)/1.e3
		# to print
		value_part   = str(res_part)+' '+units_part
		value_fields = str(res_fields)+' '+units_fields
		memory_part.set(value_part)
		memory_fields.set(value_fields)
###################
		if Interp_order_value == 'linear':
			Nshape = 6.
		elif Interp_order_value == 'quadratic':
			Nshape =19.
		FLOPns=1. # 1 double-FLoating-point instructions (operate, store or load) per nanosecond 
		# Haswell CPU architecture experimental correction factor
		FLOPns= 1.1939*FLOPns
		if Ndim_value == '1D':
			TYee        =  5.*Nx_value*Nt_value/FLOPns
			TBoris      = 23.*Npart_value*Nt_value/FLOPns
			if Interp_order_value == 'linear':
				Tdepos  = 27.*Npart_value*Nt_value/FLOPns
				Tinterp = 10.*Npart_value*Nt_value/FLOPns
			elif Interp_order_value == 'quadratic':
				Tdepos  = 40.*Npart_value*Nt_value/FLOPns
				Tinterp = 15.*Npart_value*Nt_value/FLOPns
			TEsirkepov  =  (5.*2.*Nshape+10.+5.)*Npart_value*Nt_value/FLOPns
		elif Ndim_value == '2D':
			TYee        = 35.*Nx_value*Ny_value*Nt_value/FLOPns
			TBoris      = 73.*Npart_value*Nt_value/FLOPns
			if Interp_order_value == 'linear':
				Tdepos  = 44.*Npart_value*Nt_value/FLOPns
				Tinterp = 72.*Npart_value*Nt_value/FLOPns
			elif Interp_order_value == 'quadratic':
				Tdepos  = 113.*Npart_value*Nt_value/FLOPns
				Tinterp = 135.*Npart_value*Nt_value/FLOPns
			TEsirkepov  = (10.*2.*Nshape+15.+2.*25.)*Npart_value*Nt_value/FLOPns
		elif Ndim_value == '3D':
			TYee        =  90.*Nx_value*Ny_value*Nz_value*Nt_value/FLOPns
			TBoris      = 128.*Npart_value*Nt_value/FLOPns
			if Interp_order_value == 'linear':
				Tdepos  =  60.*Npart_value*Nt_value/FLOPns
				Tinterp = 234.*Npart_value*Nt_value/FLOPns
			elif Interp_order_value == 'quadratic':
				Tdepos  = 197.*Npart_value*Nt_value/FLOPns
				Tinterp = 456.*Npart_value*Nt_value/FLOPns
			TEsirkepov  = (15.*2.*Nshape+30.+3.*125.)*Npart_value*Nt_value/FLOPns
		Ttot = TBoris + TYee + Tinterp  + Tdepos + TEsirkepov
		TBoris2     = int(1.e2*100.*TBoris/Ttot)/100
		Tdepos2     = int(1.e2*100.*Tdepos/Ttot)/100
		TEsirkepov2 = int(1.e2*100.*TEsirkepov/Ttot)/100
		TYee2       = int(1.e2*100.*TYee/Ttot)/100
		Tinterp2    = int(1.e2*100.*Tinterp/Ttot)/100
		ns          = 1./(3.6e12) # 1 ns in hours
		Ttot2       = int(1.e2*Ttot*ns)/100
		# to print
		value_comput_time     = str(Ttot2)+" CPU x hours"
		value_Boris_time      = str(TBoris2)+' %'
		value_deposit_time    = str(Tdepos2)+' %'
		value_Esirkepov_time  = str(TEsirkepov2)+' %'
		value_Yee_time        = str(TYee2)+' %'
		value_EMfields_time   = str(Tinterp2)+' %'
		comput_time.set(value_comput_time)
		Boris_time.set(value_Boris_time)
		deposit_time.set(value_deposit_time)
		Esirkepov_time.set(value_Esirkepov_time)
		Yee_time.set(value_Yee_time)
		EMfields_time.set(value_EMfields_time)
	except ValueError:
		pass

title = "Particle-In-Cell (PIC) Simulation Cost Estimator by Dr Michael J TOUATI - 06/28/2019 - mtouati@clpu.es" 

root = Tk()
root.title(title)
mainframe = ttk.Frame(root, padding="14 5 14 14")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

# Create the different widgets : 

Ndim = StringVar()
Ndim_entry = ttk.Combobox(mainframe, width=7, textvariable=Ndim,values=['1D', '2D', '3D'])
Ndim_entry.grid(column=2, row=1, sticky=(W, E))
ttk.Label(mainframe, text="Number of dimensions:").grid(column=1, row=1, sticky=W)

Interp_order = StringVar()
Interp_order_entry = ttk.Combobox(mainframe, width=7, textvariable=Interp_order,values=['linear', 'quadratic'])
Interp_order_entry.grid(column=5, row=1, sticky=(W, E))
ttk.Label(mainframe, text="Order of fields interpolation / particle shape:").grid(column=4, row=1, sticky=W)

Nspecies = StringVar()
Nspecies_entry = ttk.Entry(mainframe, width=7, textvariable=Nspecies)
Nspecies_entry.grid(column=5, row=2, sticky=(W, E))
ttk.Label(mainframe, text="Number of Species:").grid(column=4, row=2, sticky=W)

Npart = StringVar()
Npart_entry = ttk.Entry(mainframe, width=7, textvariable=Npart)
Npart_entry.grid(column=5, row=3, sticky=(W, E))
ttk.Label(mainframe, text="Total number of macro-particles (all species):").grid(column=4, row=3, sticky=W)

Nt = StringVar()
Nt_entry = ttk.Entry(mainframe, width=7, textvariable=Nt)
Nt_entry.grid(column=5, row=4, sticky=(W, E))
ttk.Label(mainframe, text="Number of time iterations:").grid(column=4, row=4, sticky=W)

Nx = StringVar()
Nx_entry = ttk.Entry(mainframe, width=7, textvariable=Nx)
Nx_entry.grid(column=2, row=2, sticky=(W, E))
ttk.Label(mainframe, text="Number of grid points along the x-axis (1D):").grid(column=1, row=2, sticky=W)

Ny = StringVar()
Ny_entry = ttk.Entry(mainframe, width=7, textvariable=Ny)
Ny_entry.grid(column=2, row=3, sticky=(W, E))
ttk.Label(mainframe, text="Number of grid points along the y-axis (2D):").grid(column=1, row=3, sticky=W)
ttk.Label(mainframe, text="(Fill 0 for 1D simulations)").grid(column=3, row=3, sticky=W)

Nz = StringVar()
Nz_entry = ttk.Entry(mainframe, width=7, textvariable=Nz)
Nz_entry.grid(column=2, row=4, sticky=(W, E))
ttk.Label(mainframe, text="Number of grid points along the z-axis (3D):").grid(column=1, row=4, sticky=W)
ttk.Label(mainframe, text="(Fill 0 for 1D or 2D simulations)").grid(column=3, row=4, sticky=W)
try:
	img = Image.open('PIC_scheme.png')
	img = ImageTk.PhotoImage(img)
	ttk.Label(mainframe, image = img).grid(row=6, column=3)
except:
    print('Image not found')
ttk.Button(mainframe, text="Roughly estimate the PIC simulation cost", command=estimate).grid(column=3, row=7, sticky=W)

ttk.Label(mainframe, text="Needed computer RAM memory:").grid(column=1, row=7, sticky=W)

ttk.Label(mainframe, text="* for macro-particles data arrays storage:").grid(column=1, row=8, sticky=W)
memory_part = StringVar()
ttk.Label(mainframe, textvariable=memory_part).grid(column=2, row=8, sticky=(W, E))

ttk.Label(mainframe, text="* for fields data arrays storage:").grid(column=1, row=9, sticky=W)
memory_fields = StringVar()
ttk.Label(mainframe, textvariable=memory_fields).grid(column=2, row=9, sticky=(W, E))

ttk.Label(mainframe, text="It has been assumed the use of:").grid(column=1, row=11, sticky=W)
ttk.Label(mainframe, text="* double-floating-point precision").grid(column=1, row=12, sticky=W)
ttk.Label(mainframe, text="* Haswell CPU architecture").grid(column=1, row=13, sticky=W)

#############

ttk.Label(mainframe, text="Needed computational time:").grid(column=4, row=7, sticky=W)

ttk.Label(mainframe, text="Total:").grid(column=4, row=8, sticky=W)
comput_time = StringVar()
ttk.Label(mainframe, textvariable=comput_time).grid(column=5, row=8, sticky=(W, E))

ttk.Label(mainframe, text="* Boris macro-particles pusher:").grid(column=4, row=9, sticky=W)
Boris_time = StringVar()
ttk.Label(mainframe, textvariable=Boris_time).grid(column=5, row=9, sticky=(W, E))

ttk.Label(mainframe, text="* Macro-particles charge deposit:").grid(column=4, row=10, sticky=W)
deposit_time = StringVar()
ttk.Label(mainframe, textvariable=deposit_time).grid(column=5, row=10, sticky=(W, E))

ttk.Label(mainframe, text="* Esirkepov charge conserving scheme:").grid(column=4, row=11, sticky=W)
Esirkepov_time = StringVar()
ttk.Label(mainframe, textvariable=Esirkepov_time).grid(column=5, row=11, sticky=(W, E))

ttk.Label(mainframe, text="* Second order Yee scheme:").grid(column=4, row=12, sticky=W)
Yee_time = StringVar()
ttk.Label(mainframe, textvariable=Yee_time).grid(column=5, row=12, sticky=(W, E))

ttk.Label(mainframe, text="* Electromagnetic fields interpolation:").grid(column=4, row=13, sticky=W)
EMfields_time = StringVar()
ttk.Label(mainframe, textvariable=EMfields_time).grid(column=5, row=13, sticky=(W, E))

############

for child in mainframe.winfo_children(): child.grid_configure(padx=14, pady=14)
Ndim_entry.focus()
Interp_order_entry.focus()
Nspecies_entry.focus()
Npart_entry.focus()
Nx_entry.focus()
Ny_entry.focus()
Nz_entry.focus()
Nt_entry.focus()
root.bind('<Return>', estimate)
root.mainloop()
