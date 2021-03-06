# a tutorial adapted from the fipy diffusion 1D example
# see: http://www.ctcms.nist.gov/fipy/examples/diffusion/index.html
# to run this file, first load JFVM and GR packages
using JFVM, GR, SpecialFunctions

function diffusion_tutorial()
    # define the domain
    L = 5  # domain length
    Nx = 100 # number of cells
    meshstruct = createMesh1D(Nx, L)
    BC = createBC(meshstruct) # all Neumann boundary condition structure
    BC.left.a[:] .= 0.0
    BC.left.b[:] .=1.0
    BC.left.c[:] .=1.0 # left boundary
    BC.right.a[:] .= 0.0
    BC.right.b[:] .=1.0
    BC.right.c[:] .=0.0 # right boundary
    x = meshstruct.cellcenters.x
    ## define the transfer coeffs
    D_val = 1.0
    alfa = createCellVariable(meshstruct, 1)
    Dave = createFaceVariable(meshstruct, D_val)
    ## define initial values
    c_old = createCellVariable(meshstruct, 0, BC) # initial values
    c = c_old
    ## loop
    dt = 0.001 # time step
    final_t = 1.5
    c_analytical=zeros(size(x))
    for t=dt:dt:final_t
        # step 1: calculate divergence term
        RHS = divergenceTerm(Dave*gradientTerm(c_old))
        c = solveExplicitPDE(c_old, dt, RHS, BC, alfa)
        # step 2: calculate the new value for internal cells
        # c_old_int=internalCells(c_old)
        # c_val_int = dt*excludeGhostRHS(meshstruct, RHS+constantSourceTerm(c_old))+c_old_int[:]
        c_analytical .= 1.0 .- erf.(x/(2*sqrt(D_val*t)))
        # GR.plot(x, [c.value[2:end-1] c_analytical])
        c_old = c
        println(t)
    end
    GR.plot(x, [c.value[2:end-1] c_analytical])
end

diffusion_tutorial()
