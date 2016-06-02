"""
# Imbibition cell, IMPES
# Buckley Leverett equation
# dependent variables: pressure and water saturation
# Prepared for educational purposes by ** AAE **
# load peteng relperm files before running this function
## define the geometry
"""
function imb_impes()
Nx = 20 # number of cells in x direction
Ny = 50 # number of cells in y direction
W = 0.02 # [m] length of the domain in x direction
H = 0.07 # [m] length of the domain in y direction
# m = createMesh1D(Nx, W)
m = createMesh2D(Nx, Ny, W, H) # creates a 2D mesh
## define the physical parametrs
# all the Corey-type relperm parametrs are defined for the oil-wet and water-wet
# cases
k0 = 0.01e-12 # [m^2] average reservoir permeability
phi0 = 0.45 # average porosity
mu_oil = 2e-3 # [Pa.s] oil viscosity
mu_water = 1e-3 # [Pa.s] water viscosity
krw0_ww = 0.3
krw0_ow = 1.0
kro0_ww = 0.6
kro0_ow = 0.76
nw_ww = 2.4
nw_ow= 2.4
no_ww = 2.0
no_ow= 2.0
sor_ww=0.1
sor_ow=0.12
swc_ww=0.09
swc_ow=0.09
teta_ow=deg2rad(130)
teta_ww=deg2rad(20)
gama_ow=0.03 # N/m
labda=2.4 # for Ekofisk chalk
r_ave=sqrt(k0/phi0) #  meter average pore diameter
pce0=2*gama_ow*cos(teta_ww)/r_ave # Pa capillary entry pressure

SF=createFaceVariable(m, 0.0) # 1 is water wet, 0 is oil wet
krw0=krw0_ww*SF+krw0_ow*(1-SF)
kro0=kro0_ww*SF+kro0_ow*(1-SF)
sor=sor_ww*SF+sor_ow*(1-SF)
swc=swc_ww*SF+swc_ow*(1-SF)
no= no_ww*SF+no_ow*(1-SF)
nw= nw_ww*SF+nw_ow*(1-SF)
teta=teta_ww*SF+teta_ow*(1-SF)
pce=createFaceVariable(m, pce0)

p0 = 100e5 # [bar] pressure
pin = 150e5 # [bar] injection pressure at the left boundary
u_in= 1.0/(24*3600) # [m/s] equal to 1 m/day
sw0 = swc_ww+0.01
# sw0(10:end-10, 10:end-10)=swc+0.2
# sw0 = swc+0.1 # initial water saturation
sw_in = 1

eps1=1e-6
clx=0.05
cly=0.05
V_dp=0.01 # Dykstra-Parsons coef.
perm_val= k0 #field2d(Nx,Ny,k0,V_dp,clx,cly)

k=createCellVariable(m, perm_val)
phi=createCellVariable(m, phi0)

lw = geometricMean(k)/mu_water
lo = geometricMean(k)/mu_oil

## Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
BCp = createBC(m) # Neumann BC for pressure
BCs = createBC(m) # Neumann BC for saturation
# left boundary pressure gradient
# BCp.left.a[:]=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)) BCp.left.b[:]=0 BCp.left.c[:]=-u_in
# change the right boandary to constant pressure (Dirichlet)
# BCp.left.a[:]=0 BCp.left.b[:]=1 BCp.left.c[:]=p0
BCp.right.a[:]=0.0
BCp.right.b[:]=1.0
BCp.right.c[:]=p0
BCp.top.a[:]=0.0
BCp.top.b[:]=1.0
BCp.top.c[:]=p0
BCp.bottom.a[:]=0.0
BCp.bottom.b[:]=1.0
BCp.bottom.c[:]=p0
# change the left boundary to constant saturation (Dirichlet)
# BCs.left.a[:]=0 BCs.left.b[:]=1 BCs.left.c[:]=1.0-sor
BCs.right.a[:]=0.0
BCs.right.b[:]=1.0
BCs.right.c[:]=1.0
BCs.top.a[:]=0.0
BCs.top.b[:]=1.0
BCs.top.c[:]=1.0
BCs.bottom.a[:]=0.0
BCs.bottom.b[:]=1.0
BCs.bottom.c[:]=1.0
## define the time step and solver properties
# dt = 1000 # [s] time step
# dt=(W/Nx)/u_in/20 # [s]
dt=1
t_end = 10*3600*24 # [s] final time
eps_p = 1e-7 # pressure accuracy
eps_sw = 1e-7 # saturation accuracy
## define the variables
sw_old = createCellVariable(m, sw0, BCs)
p_old = createCellVariable(m, p0, BCp)
sw = copyCell(sw_old)
oil_init=domainInt(1-sw_old)
p = copyCell(p_old)
uw = -gradientTerm(p_old) # an estimation of the water velocity
## start the main loop
# generate intial pressure profile (necessary to initialize the fully
# implicit solver)
rec_fact=zeros(1)
t_day=zeros(1)
t = 0.0
dt0=dt
dsw_alwd= 0.001
dp_alwd= 100.0 # Pa
while (t<t_end)
# for i=1:5
    error_p = 1e5
    error_sw = 1e5
    # Implicit loop
#     while ((error_p>eps_p) || (error_sw>eps_sw))
    while(true)
        # calculate parameters
        pgrad = gradientTerm(p)
#         pcgrad=gradientTerm(pc(sw))
        sw_face = upwindMean(sw, -pgrad) # average value of water saturation
        sw_grad=gradientTerm(sw)
        sw_ave=arithmeticMean(sw)
        pcgrad=faceEval(dpc_imb, sw_ave, pce, swc, sor, teta)
        # solve for pressure at known Sw
        labdao = lo.*faceEval(kro, sw_face, kro0, sor, swc, no)
        labdaw = lw.*faceEval(krw, sw_face, krw0, sor, swc, nw)
        labda = labdao+labdaw
        # compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda)
        RHSpc1=divergenceTerm(labdao.*pcgrad)
        Mbcp, RHSbcp = boundaryConditionTerm(BCp)
        RHS1 = RHSpc1+RHSbcp # with capillary
        p_new=solvePDE(m, Mdiffp1+Mbcp, RHS1)

        # solve for Sw
        pgrad = gradientTerm(p_new)
        uw=-labdaw.*pgrad
        Mbcsw, RHSbcsw = boundaryConditionTerm(BCs)
        RHS_sw=-divergenceTerm(uw)
        sw_new=solveExplicitPDE(sw_old, dt, RHS_sw, BCs, phi)

        error_p = maximum(abs((internalCells(p_new)-internalCells(p))./internalCells(p_new)))
        error_sw = maximum(abs(internalCells(sw_new)-internalCells(sw)))
        dt_new=dt*min(dp_alwd/error_p, dsw_alwd/error_sw)
        # print("sw_error = $error_sw \n")
        # print("new time step = $dt \n")
        # print(internalCells(sw_new))
        # print(internalCells(sw_old))
        # sleep(1.0)
        # assign new values of p and sw
        if error_sw>dsw_alwd
            dt=dt*(dsw_alwd/error_sw)
            # print("new time step = $dt \n")
        else
            t=t+dt
            p = copyCell(p_new)
            sw = copyCell(sw_new)
            p_old = copyCell(p)
            sw_old = copyCell(sw)
            dt=min(dt*(dsw_alwd/error_sw), 10*dt)
            break
        end
    end

    rec_fact=push!(rec_fact, (oil_init-domainInt(1-sw))/oil_init)
    t_day=push!(t_day, t)
    # print(t)
    GR.imshow(sw.value[2:end-1,2:end-1])
    # visualizeCells(1-sw)
    # GR.plot(t_day/3600/24, rec_fact)
    # xlabel('time [day]')
    # ylabel('recovery factor')
    # title([num2str(t/3600/24) ' day']) drawnow
end
end # imb_impes
