
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
from yt.mods import *
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import ImageGrid


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# nyx resolutions
dims = np.array([32,64,128,256,512])#np.array([32,64,128,256,512])
# plotfile numbers
pfs = np.array([[81,87,92,96,101,104,107],
                [81,87,92,96,101,105,111],
                [81,88,95,103,116,127,142],
                [86,100,116,134,164,189,225],
                [113,142,178,220,289,347,431]])
# corresponding z's
z = [20.,15.,12.,10.,8.,7.,6.]

# Density and Temperature slices, fixed redshift and varying resolution

eta2 = "0.001"

fig = plt.figure(figsize=(3.0*len(dims),9))

iz = -2

# Load in density and temperature fields
dens = [None for ii in range(len(dims))]
temp = [None for ii in range(len(dims))]
mach = [None for ii in range(len(dims))]
for ii in range(len(dims)):
    pf = load('{:d}_e{:s}/plt{:05d}'.format(dims[ii],eta2,pfs[ii][iz]))
    dens[ii] = pf.covering_grid(0,left_edge = [0.0,0.0,0.0],dims=[dims[ii],dims[ii],dims[ii]],fields=["density"])['density']
    dens[ii] = dens[ii]/np.mean(dens[ii])
    temp[ii] = pf.covering_grid(0,left_edge = [0.0,0.0,0.0],dims=[dims[ii],dims[ii],dims[ii]],fields=["Temp"])['Temp']
    mach[ii] = pf.covering_grid(0,left_edge = [0.0,0.0,0.0],dims=[dims[ii],dims[ii],dims[ii]],fields=["MachNumber"])['MachNumber']
#    dx = pf.all_data()['dx'][0]
#    x = np.floor((pf.all_data()['x']/dx).value).astype(int)
#    y = np.floor((pf.all_data()['y']/dx).value).astype(int)
#    z = np.floor((pf.all_data()['z']/dx).value).astype(int)
#    linear_index = z+dims[ii]*y+dims[ii]*dims[ii]*x
#    dens[ii] = pf.all_data()['density'][linear_index].reshape((dims[ii],dims[ii],dims[ii]))
#    dens[ii] = dens[ii]/np.mean(dens[ii])
#    temp[ii] = pf.all_data()['Temp'][linear_index].reshape((dims[ii],dims[ii],dims[ii]))

dens_range = [-1.0,1.0]#[-1,1]
temp_range = [-2,3.0]#[-2.5,4]
mach_range = [0,2.8]

grid = ImageGrid(fig,111,nrows_ncols=(3,len(dims)),axes_pad=(0.001,0.08),cbar_mode='edge',cbar_pad=0.06,cbar_location='right')

for ii in range(len(dims)):
    slice = dims[ii]//2-1
    imd = grid[ii].imshow(np.log10(dens[ii][:,slice,:]),vmin=dens_range[0],vmax=dens_range[1],
                          cmap='BrBG_r',extent=[0,1,0,1],aspect='equal')
    grid[ii].set_title(r'{:d}$^3$, $z=$ {:.1f}, $\eta_2=$ {:s}'.format(dims[ii],z[iz],eta2),fontsize=15)
    grid[ii].tick_params(which='both',direction='in',top=True,right=True,
                         labelbottom=False,labelleft=False)
    grid[ii].minorticks_on()
    imt = grid[ii+len(dims)].imshow(np.log10(temp[ii][:,slice,:]),vmin=temp_range[0],vmax=temp_range[1],
                            cmap='RdBu_r',extent=[0,1,0,1],aspect='equal')
    grid[ii+len(dims)].tick_params(which='both',direction='in',top=True,right=True,
                         labelbottom=False,labelleft=False)
    grid[ii+len(dims)].minorticks_on()
    imm = grid[ii+2*len(dims)].imshow(np.log10(mach[ii][:,slice,:]),vmin=mach_range[0],vmax=mach_range[1],
                                    cmap='Spectral_r',extent=[0,1,0,1],aspect='equal')
    grid[ii+2*len(dims)].tick_params(which='both',direction='in',top=True,right=True,
                                   labelbottom=False,labelleft=False)
    grid[ii+2*len(dims)].minorticks_on()


cb = grid.cbar_axes[0].colorbar(imd)
cb.set_label_text(r'log $\Delta_b$',fontsize=14)
grid.cbar_axes[0].minorticks_on()
cb2 = grid.cbar_axes[1].colorbar(imt)
cb2.set_label_text(r'log $T$',fontsize=14)
grid.cbar_axes[1].minorticks_on()
cb3 = grid.cbar_axes[2].colorbar(imm)
cb3.set_label_text(r'log $\mathcal{M}$',fontsize=14)
grid.cbar_axes[2].minorticks_on()

grid[0].set_ylabel('Baryon density',fontsize=15)
grid[len(dims)].set_ylabel('Temperature',fontsize=15)
grid[len(dims)*2].set_ylabel('Mach number',fontsize=15)

plt.savefig('dt_slices_z{:.0f}_e{:s}.pdf'.format(z[iz],eta2),dpi=350)

# T-rho histograms, fixed redshift varying resolution

ranges = [[[-0.45,0.75],[0.0,2.5]],
          [[-0.5,0.75],[-0.5,3.0]],
          [[-0.75,1.0],[-1.0,3.0]],
          [[-1.0,2.0],[-1.0,4.0]],
          [[-0.9,2.0],[-2.0,4.0]],
          [[-0.9,2.0],[-2.5,4.0]],
          [[-0.9,2.0],[-2.5,4.0]]]

fig,ax = plt.subplots(1,len(dims),figsize=(3*len(dims),4))

rf = np.loadtxt('RECFAST')
rf_interp = interp1d(rf[:,0],np.log10(rf[:,2]))
iz = -2
d = np.arange(ranges[iz][0][0]+0.1*(ranges[iz][0][1]-ranges[iz][0][0]),
              ranges[iz][0][1]-0.1*(ranges[iz][0][1]-ranges[iz][0][0]),0.01)
t = np.log10(10**rf_interp(z[iz])*(10**d)**(2.0/3.0))
#t2 = np.log10(10**rf_interp(z[iz])*(2**(2.0/3.0))*((10**d)/2.0)**(3))
for ii in range(len(dims)):
    #pf = load('{:d}/plt{:05d}'.format(dims[ii],pfs[ii][iz]))
    
    #dens = pf.all_data()['density'].reshape((dims[ii],dims[ii],dims[ii]))
    #dens /= np.mean(dens)
    #temp = pf.all_data()['Temp'].reshape((dims[ii],dims[ii],dims[ii]))

    ax[ii].hist2d(np.log10(dens[ii]).flatten(),np.log10(temp[ii]).flatten(),bins=150,range=ranges[iz],
                  normed=True,cmap='Spectral_r',cmin=0.001)
    ax[ii].set_xlabel('log $\Delta_b$',fontsize=14)
    ax[ii].set_title('{:d}$^3$, $z=${:.1f}, $\eta_2=$ {:s}'.format(dims[ii],z[iz],eta2),fontsize=14)
    ax[ii].plot(d,t,c='k',lw=2.0,linestyle='dashed')
#ax[ii].plot(d,t2,c='k',lw=2.0,linestyle='dotted')

ax[0].set_ylabel('log $T$',fontsize=14)

for x in ax.flatten():
    x.tick_params(which='both',direction='in',top=True,right=True,
                  labelleft=False)
    x.minorticks_on()

ax[0].tick_params(labelleft=True)

plt.tight_layout()
plt.subplots_adjust(wspace=0.001)

plt.savefig('dt_2dhists_z{:.0f}_e{:s}.pdf'.format(z[iz],eta2),dpi=200)

exit(0)

# Temperature histograms

ranges = [[-0.5,2.5],[-1.5,3.0],[-2.5,3.5],[-2.5,3.5],[-2.5,3.5]]

fig = plt.figure(figsize=(6,5))

iz = -1
for ii in range(len(dims)):
    pf = load('{:d}/plt{:05d}'.format(dims[ii],pfs[ii][iz]))
    temp = pf.all_data()['Temp']
    plt.hist(np.log10(temp),bins=50,range=ranges[iz],density=True,label='Nyx {:d}'.format(dims[ii]),alpha=0.35)
    plt.axvline(rf_interp(z[iz]),c='k',linestyle='dashed',lw=2.0)

plt.tick_params(which='both',direction='in',right=True,top=True,labelsize=13)
plt.minorticks_on()
#plt.ylim(bottom=1e-2)
#plt.yscale('log')
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

# Temperature evolution curves, from the runlog

fig = plt.figure(figsize=(6,5))

rf = np.loadtxt('RECFAST')
cmb = 2.725*(1+rf[:,0])

nt = 8

for ii in range(len(dims)-1):
    rl = np.loadtxt('outputs_from_igm/{:d}/runlog'.format(dims[ii]))
    plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**ii,linestyle='dotted')#,label=r'1 Mpc/$h$ box, {:.1f} kpc res.'.format(1000/0.675/dims[ii]))

rl = np.loadtxt('outputs_from_igm/256/runlog'.format(dims[ii]))
plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**3,linestyle='dotted',label=r'1 Mpc/$h$ box, original'.format(1000/0.675/512))

for ii in range(len(dims)-1):
    rl = np.loadtxt('outputs_from_igm/e0.03/{:d}/runlog'.format(dims[ii]))
    plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**ii)#,label=r'1 Mpc/$h$ box, {:.1f} kpc res.'.format(1000/0.675/dims[ii]))

rl = np.loadtxt('outputs_from_igm/e0.03/256/runlog'.format(dims[ii]))
plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**3,label=r'1 Mpc/$h$ box, $\eta_2=0.1$'.format(1000/0.675/512))

#rl = np.loadtxt('outputs_from_igm/512/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**4,label=r'1 Mpc/$h$ box, {:.1f} kpc res.'.format(1000/0.675/512))

#rl = np.loadtxt('outputs_from_igm/32_nominxe/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.0,linestyle='dotted')
#
#rl = np.loadtxt('outputs_from_igm/64_nominxe/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3,linestyle='dotted',label=r'No minimum $x_e$')
#
#rl = np.loadtxt('outputs_from_igm/128_nominxe/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3*1.3,linestyle='dotted')
#
#rl = np.loadtxt('outputs_from_igm/256_nominxe/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3*1.3*1.3,linestyle='dotted')
#
#rl = np.loadtxt('outputs_from_igm/512_nominxe/runlog'.format(dims[ii]))
#plt.plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3*1.3*1.3*1.3,linestyle='dotted')

plt.plot(1+rf[:,0],rf[:,2],c='darkorange',linestyle='dashed',lw=2.0,label='RECFAST')
plt.plot(1+rf[:,0],cmb,c='dodgerblue',linestyle='dashed',lw=2.0,label='CMB')

plt.xscale('log')
plt.yscale('log')

plt.gca().set_yticklabels([0,0,1,10,100,1000])
plt.gca().set_xticklabels([0,0,10])
plt.gca().set_xticklabels([0,0,0,0,0,0,0,0,0,0,0,0,0,"",8,9,20,30,40],minor=True)

plt.tick_params(which='both',direction='in',right=True,top=True,labelsize=13)
plt.xlabel('1+$z$',fontsize=20)
plt.ylabel(r'$T(\langle \rho \rangle)$ [K]',fontsize=20)
plt.xlim(7,51)
plt.ylim(0.5,3000)

plt.legend(fontsize=12)

plt.tight_layout()
plt.show()





etas = ["0.3","0.1","0.01","0.001"]

fig,ax = plt.subplots(1,len(etas),figsize=(4*len(etas),5))

rf = np.loadtxt('RECFAST')
cmb = 2.725*(1+rf[:,0])

nt = 9

for ie in range(len(etas)):
    for ii in range(len(dims)):
        rl = np.loadtxt('outputs_from_igm/runlogs/runlog_{:d}_e{:s}'.format(dims[ii],etas[ie]))
        if ie == 0:
            ax[ie].plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**ii,label=r'{:.1f} kpc res.'.format(1000/0.675/dims[ii]))
        else:
            ax[ie].plot(1+rl[:,3],rl[:,nt],c='k',lw=1.3**ii)
    if ie == 0:
        ax[ie].plot(1+rf[:,0],rf[:,2],c='darkorange',linestyle='dashed',lw=2.0,label='RECFAST')
        ax[ie].plot(1+rf[:,0],cmb,c='dodgerblue',linestyle='dashed',lw=2.0,label='CMB')
    else:
        ax[ie].plot(1+rf[:,0],rf[:,2],c='darkorange',linestyle='dashed',lw=2.0)
        ax[ie].plot(1+rf[:,0],cmb,c='dodgerblue',linestyle='dashed',lw=2.0)
    ax[ie].set_title('$\eta_2=$ {:s}'.format(etas[ie]),fontsize=20)

for x in ax:
    x.set_xscale('log')
    x.set_yscale('log')

    x.set_yticklabels([0,0,1,10,100])
    x.set_xticklabels([0,0,10])
    x.set_xticklabels([0,0,0,0,0,0,0,0,0,0,0,0,0,"",8,9,20,30,40],minor=True)

    x.tick_params(which='both',direction='in',labelleft=False,right=True,top=True,labelsize=14)
    x.set_xlabel(r'$1+z$',fontsize=20)
    
    x.set_xlim(7,51)
    x.set_ylim(0.2,200)

if nt == 8:
    ax[0].set_ylabel(r'$T(\langle \rho \rangle)$ [K]',fontsize=20)
elif nt == 9:
    ax[0].set_ylabel(r'$T_{\rm 21 cm}$ [K]',fontsize=20)
ax[0].tick_params(which='both',labelleft=True)

ax[0].legend(fontsize=11)

plt.tight_layout()
plt.subplots_adjust(wspace=0.001)
plt.savefig('T_21cm.pdf')
plt.show()





# T-rho evolution

#ranges = [[[-0.45,0.75],[0.0,2.5]],
#          [[-0.5,0.75],[-0.5,3.0]],
#          [[-0.75,1.0],[-1.0,3.0]],
#          [[-1.0,2.0],[-2.0,4.0]],
#          [[-1.0,2.0],[-2.0,4.0]],
#          [[-1.0,2.0],[-2.5,4.0]]]
#
#fig,ax = plt.subplots(1,len(dims),figsize=(3.5*len(dims),4))
#
#for ii in range(len(dims)):
#    for iz in range(len(z)):
#    pf = load('{:d}/plt{:05d}'.format(dims[ii],pfs[ii][iz]))
#
#    dens = pf.all_data()['density'].reshape((dims[ii],dims[ii],dims[ii]))
#    dens /= np.mean(dens)
#    temp = pf.all_data()['Temp'].reshape((dims[ii],dims[ii],dims[ii]))
#
#    ax[ii].hist2d(np.log10(dens).flatten(),np.log10(temp).flatten(),bins=200,range=ranges[iz],
#                  normed=True,cmap='Spectral_r',cmin=0.001)
#                  ax[ii].set_xlabel('log10 Delta',fontsize=14)
#                  ax[ii].set_title('{:d}, z={:.1f}'.format(dims[ii],z[iz]),fontsize=14)
#                  ax[ii].plot(d,t,c='k',lw=2.0,linestyle='dashed')
##ax[ii].plot(d,t2,c='k',lw=2.0,linestyle='dotted')
#
#ax[0].set_ylabel('log10 Temperature',fontsize=14)
#
#for x in ax.flatten():
#    x.tick_params(which='both',direction='in',top=True,right=True,
#                  labelleft=False)
#    x.minorticks_on()
#
#ax[0].tick_params(labelleft=True)
#
#plt.tight_layout()
#plt.subplots_adjust(wspace=0.001)
#
#plt.show()
