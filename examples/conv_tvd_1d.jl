using JFVM, Plots


function compare_upwind_tvd()
  Nx=200
  Lx=1.0
  u=1.0
  dt=Lx/u/Nx/10
  t_end=0.5*Lx/u
  m=createMesh1D(Nx, Lx)
  u_face=createFaceVariable(m, [u])
  c_init=createCellVariable(m, 0.0)
  c_init_up=createCellVariable(m, 0.0)
  c=createCellVariable(m, 0.0)
  c_upwind=createCellVariable(m, 0.0)

  BC=createBC(m)
  BC.left.a .= 0.0
  BC.left.b .= 1.0
  BC.left.c .= 1.0
  Mbc, RHSbc=boundaryConditionTerm(BC)

  Mconv_up=convectionUpwindTerm(u_face)

  flux_limiter=fluxLimiter("SUPERBEE")
  for t in dt:dt:t_end
    Mt, RHSt=transientTerm(c_init, dt)
    Mt_up, RHSt_up=transientTerm(c_init_up, dt)
    c_upwind=solveLinearPDE(m, Mt_up+Mconv_up+Mbc, RHSt_up+RHSbc)
    for j in 1:10
      RHSconv=convectionTvdRHS(u_face, c, flux_limiter)
      c=solveLinearPDE(m, Mt+Mconv_up+Mbc, RHSt+RHSbc+RHSconv)
    end
    c_init=copyCell(c)
    c_init_up=copyCell(c_upwind)
  end
  Plots.plot(c.value[2:end-1])
  Plots.plot!(c_upwind.value[2:end-1])
end

compare_upwind_tvd()