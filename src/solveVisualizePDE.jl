# ===============================
# Written by AAE
# TU Delft, Winter 2014
# simulkade.com
# Last edited: 29 December, 2014
# ===============================

# =========================================================================
# Changes
# 2014-12-31: - Added 2D radial surface plot
# =========================================================================

# =============================== SOLVERS ===================================
function solveLinearPDE(m::MeshStructure, M::SparseMatrixCSC{Float64, Int64}, RHS::Array{Float64,1})
N=m.numberofcells
x=M\RHS
phi = CellValue(m, reshape(full(x), tuple(N+2...)))
phi
end



# =========================== Visualization =================================
function visualizeCells(phi::CellValue)
d=phi.domain.dimension
if d==1 || d==1.5
  x = [phi.domain.facecenters.x[1]; phi.domain.cellcenters.x; phi.domain.facecenters.x[end]]
  phi = [0.5*(phi.value[1]+phi.value[2]); phi.value[2:end-1]; 0.5*(phi.value[end-1]+phi.value[end])]
  plot(x, phi)
elseif d==2 || d==2.5
  x = [phi.domain.facecenters.x[1]; phi.domain.cellcenters.x; phi.domain.facecenters.x[end]]
  y = [phi.domain.facecenters.y[1]; phi.domain.cellcenters.y; phi.domain.facecenters.y[end]]
  phi0 = Base.copy(phi.value)
  phi0[:,1] = 0.5*(phi0[:,1]+phi0[:,2])
  phi0[1,:] = 0.5*(phi0[1,:]+phi0[2,:])
  phi0[:,end] = 0.5*(phi0[:,end]+phi0[:,end-1])
  phi0[end,:] = 0.5*(phi0[end,:]+phi0[end-1,:])
  phi0[1,1] = phi0[1,2] 
  phi0[1,end] = phi0[1,end-1]
  phi0[end,1] = phi0[end,2]
  phi0[end,end] = phi0[end,end-1]
  pcolor(x, y, phi0')
elseif d==2.8
  x = [phi.domain.facecenters.x[1]; phi.domain.cellcenters.x; phi.domain.facecenters.x[end]]
  y = [phi.domain.facecenters.y[1]; phi.domain.cellcenters.y; phi.domain.facecenters.y[end]]
  phi0 = Base.copy(phi.value)
  phi0[:,1] = 0.5*(phi0[:,1]+phi0[:,2])
  phi0[1,:] = 0.5*(phi0[1,:]+phi0[2,:])
  phi0[:,end] = 0.5*(phi0[:,end]+phi0[:,end-1])
  phi0[end,:] = 0.5*(phi0[end,:]+phi0[end-1,:])
  phi0[1,1] = phi0[1,2] 
  phi0[1,end] = phi0[1,end-1]
  phi0[end,1] = phi0[end,2]
  phi0[end,end] = phi0[end,end-1]
  subplot(111, polar="true")
  pcolor(y, x, phi0)  
end
end  