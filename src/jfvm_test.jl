function JFVM_test()
# JFVM test script
# This script is supposed to call all the functions of the JFVM package
## Part I: creating an array of different mesh types:
# domain size
Lx= 1.0
Ly= 2*pi
Lz= 2.0
Nx=5
Ny=7
Nz=9
X=[0.01, 0.1, 0.3, 0.5, 0.55, 1.0]
Y= [0.0, 0.1, 1.0, 1.5, 2.9, 3.0, pi, 2*pi]
Z= [0.0, 0.01, 0.1, 0.5, 0.7, 0.95, 1.0, 1.25, 1.39, 2.0]
N_mesh=7
# create nonuniform mesh
mesh_nonuniform= Array{Any}(undef, N_mesh)
mesh_nonuniform[1]=createMesh1D(X)
mesh_nonuniform[2]=createMesh2D(X, Y)
mesh_nonuniform[3]=createMesh3D(X, Y, Z)
mesh_nonuniform[4]=createMeshCylindrical1D(X)
mesh_nonuniform[5]=createMeshCylindrical2D(X, Y)
mesh_nonuniform[6]=createMeshCylindrical3D(X, Y, Z)
mesh_nonuniform[7]=createMeshRadial2D(X, Y)
println("Non-uniform mesh created successfully!")
## Part II: create cell and face variables
c_val= 1.0
# nonuniform
c_n=Array{Any}(undef, N_mesh)
c_r=Array{Any}(undef, N_mesh)
D_n=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    c_n[i]= createCellVariable(mesh_nonuniform[i], c_val)
end
for i=1:N_mesh
    D_n[i]= createCellVariable(mesh_nonuniform[i], c_val)
end
println("Cells of fixed values over nonuniform mesh created successfully!")
for i=1:N_mesh
    c_rand = rand(Float64, tuple(mesh_nonuniform[i].dims...))
    c_r[i]= createCellVariable(mesh_nonuniform[i], c_rand)
end
println("Cells of random values over nonuniform mesh created successfully!")
## Part III: create face variables
f_val= 0.5
# nonuniform
f_n=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    f_n[i]= createFaceVariable(mesh_nonuniform[i], f_val)
end
println("Face variable over nonuniform mesh created successfully!")
## Part IV: Test boundary conditions
BC_n=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    BC_n[i]=createBC(mesh_nonuniform[i])
    BC_n[i].left.a[:].=0.0
    BC_n[i].left.b[:].=1.0
    BC_n[i].left.c[:].=1.0
    BC_n[i].right.a[:].=0.0
    BC_n[i].right.b[:].=1.0
    BC_n[i].right.c[:].=0.0
end
println("Boundary condition over a nonuniform mesh created successfully!")
## Part V: solve a steady-state diffusion equation
c_dif=Array{Any}(undef, N_mesh)
M_bc=Array{Any}(undef, N_mesh)
M_dif=Array{Any}(undef, N_mesh)
RHS_bc=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    M_bc[i], RHS_bc[i]= boundaryConditionTerm(BC_n[i])
    M_dif[i]=diffusionTerm(f_n[i])
    c_dif[i]=solveLinearPDE(mesh_nonuniform[i], M_dif[i]+M_bc[i], RHS_bc[i])
end
# visualize
# figure(1)
# for i=1:N_mesh
#     subplot(3, 3, i)
#     visualizeCells(c_dif[i])
# end
# println("Diffusion equation solved and visualized successfully")
## Part VI: solve convection diffucion equation
# nonuniform
c_conv=Array{Any}(undef, N_mesh)
M_bc=Array{Any}(undef, N_mesh)
M_dif=Array{Any}(undef, N_mesh)
M_conv=Array{Any}(undef, N_mesh)
RHS_bc=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    M_bc[i], RHS_bc[i]= boundaryConditionTerm(BC_n[i])
    M_dif[i]=diffusionTerm(f_n[i])
    M_conv[i]=convectionTerm(0.1*f_n[i])
    c_conv[i]=solveLinearPDE(mesh_nonuniform[i], M_conv[i]-M_dif[i]+M_bc[i], RHS_bc[i])
end
# visualize
# figure(2)
# for i=1:N_mesh
#     subplot(3, 3, i)
#     visualizeCells(c_conv[i])
# end
# println("Convection-Diffusion equation solved and visualized successfully")
## Part VII: test the calculus fanctions
grad_c=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    grad_c[i]=gradientTerm(c_dif[i])
end
div_c=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    div_c[i]=divergenceTerm(grad_c[i])
end
println("Gradient and divergence functions work fine!")
## Solve a dynamic equation
dt=0.1
FL1=fluxLimiter("SUPERBEE")
c_old=c_n
c_trans=Array{Any}(undef, N_mesh)
M_bc=Array{Any}(undef, N_mesh)
M_dif=Array{Any}(undef, N_mesh)
M_conv=Array{Any}(undef, N_mesh)
RHS_bc=Array{Any}(undef, N_mesh)
M_ls=Array{Any}(undef, N_mesh)
RHS_s=Array{Any}(undef, N_mesh)
RHS_tvd=Array{Any}(undef, N_mesh)
for i=1:N_mesh
    M_bc[i], RHS_bc[i]= boundaryConditionTerm(BC_n[i])
    M_dif[i]=diffusionTerm(0.1*f_n[i])
    M_conv[i]=convectionUpwindTerm(0.01*f_n[i])
    RHS_tvd[i]=convectionTvdRHS(0.01*f_n[i], c_old[i], FL1) #only called, not used
    M_ls[i]=linearSourceTerm(0.1*c_n[i])
    RHS_s[i]=constantSourceTerm(0.2*c_n[i])
end

for i=1:N_mesh
    for j=1:10
        (M_t, RHS_t)=transientTerm(c_old[i], dt)
        c_trans[i]=solveLinearPDE(mesh_nonuniform[i],
            M_t+M_ls[i]+M_conv[i]-M_dif[i]+M_bc[i], RHS_t+RHS_s[i]+RHS_bc[i])
        c_old[i]=copyCell(c_trans[i])
    end
end
# visualize
# figure(3)
#for i=1:N_mesh
#     visualizeCells(c_trans[i])
#     pause(1.5)
#end
println("Transient convection-diffucion-reaction solved successfully!")
## Part VIII: test the utilities
# only test the averaging, don"t save the result
for i=1:N_mesh
    linearMean(c_trans[i])
    arithmeticMean(c_trans[i])
    geometricMean(D_n[i])
    harmonicMean(D_n[i])
    upwindMean(c_trans[i], f_n[i])
    tvdMean(c_trans[i], f_n[i], FL1)
end
println("Averaging functions run smoothly!")
## Part IX: test the classes and operators

end # end function
