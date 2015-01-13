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
N=m.dims
x=lufact(M)\RHS # until the problem is solved with Julia "\" solver
phi = CellValue(m, reshape(full(x), tuple(N+2...)))
phi
end

function solveMUMPSLinearPDE(m::MeshStructure, M::SparseMatrixCSC{Float64, Int64}, RHS::Array{Float64,1})
N=m.dims
x=solveMUMPS(M,RHS) # until the problem is solved with Julia "\" solver
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
elseif d==3
  x=phi.domain.cellcenters.x
  y=phi.domain.cellcenters.y
  z=phi.domain.cellcenters.z
  Nx = phi.domain.dims[1]
  Ny = phi.domain.dims[2]
  Nz = phi.domain.dims[3]
  vmin = minimum(phi.value[2:end-1,2:end-1,2:end-1])
  vmax = maximum(phi.value[2:end-1,2:end-1,2:end-1])
  # 6 surfaces
  # surfaces 1,2 (x=x[1], x=x[end])
  Y=repmat(y,1,Nz)
  Z=repmat(z,1,Ny)
  mayavis.mesh(x[1]*ones(Ny,Nz),Y,Z',scalars=squeeze(phi.value[2,2:end-1,2:end-1],1), vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(x[end]*ones(Ny,Nz),Y,Z',scalars=squeeze(phi.value[end-1,2:end-1,2:end-1],1), vmin=vmin, vmax=vmax, opacity=0.8)

  # surfaces 3,4 (y=y[1], y=y[end]
  X = repmat(x,1,Nz)
  Z = repmat(z,1,Nx)
  mayavis.mesh(X,y[1]*ones(Nx,Nz),Z',scalars=squeeze(phi.value[2:end-1,2,2:end-1],2), vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(X,y[end]*ones(Nx,Nz),Z',scalars=squeeze(phi.value[2:end-1,end-1,2:end-1],2), vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.axes()

  # surfaces 5,6 (z=z[1], z=z[end]
  X = repmat(x,1,Ny)
  Y = repmat(y,1,Nx)
  mayavis.mesh(X,Y',z[1]*ones(Nx,Ny),scalars=phi.value[2:end-1,2:end-1,2], vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(X,Y',z[end]*ones(Nx,Ny),scalars=phi.value[2:end-1,2:end-1,end-1], vmin=vmin, vmax=vmax, opacity=0.8)

  mayavis.colorbar()
  mayavis.axes()
  return mayavis.screenshot()
  
elseif d==3.2
  Nx = phi.domain.dims[1]
  Ny = phi.domain.dims[2]
  Nz = phi.domain.dims[3]
  r=phi.domain.cellcenters.x
  theta = Array(Float64,1,Ny)
  theta[:]=phi.domain.cellcenters.y
  z=Array(Float64,1,1,Nz)
  z[:]=phi.domain.cellcenters.z
  a=ones(Nx,Ny,Nz)
  R=r.*a
  TH = theta.*a
  Z = z.*a
  
  X=R.*cos(TH)
  Y=R.*sin(TH)
  
  vmin = minimum(phi.value[2:end-1,2:end-1,2:end-1])
  vmax = maximum(phi.value[2:end-1,2:end-1,2:end-1])
  # 6 surfaces
  # surfaces 1,2 (x=x[1], x=x[end])
  mayavis.mesh(squeeze(X[floor(Nx/2.0),:,:],1),squeeze(Y[floor(Nx/2.0),:,:],1),squeeze(Z[floor(Nx/2.0),:,:],1), 
    scalars=squeeze(phi.value[floor(Nx/2.0)+1,2:end-1,2:end-1],1), vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(squeeze(X[Nx,:,:],1),squeeze(Y[Nx,:,:],1),squeeze(Z[Nx,:,:],1), 
    scalars=squeeze(phi.value[Nx+1,2:end-1,2:end-1],1), vmin=vmin, vmax=vmax, opacity=0.8)
    
  # surfaces 3,4 (y=y[1], y=y[end]
  mayavis.mesh(squeeze(X[:,floor(Ny/2.0),:],2),squeeze(Y[:,floor(Ny/2.0),:],2),squeeze(Z[:,floor(Ny/2.0),:],2), 
    scalars=squeeze(phi.value[2:end-1,floor(Ny/2.0)+1,2:end-1],2), vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(squeeze(X[:,Ny,:],2),squeeze(Y[:,Ny,:],2),squeeze(Z[:,Ny,:],2), 
    scalars=squeeze(phi.value[2:end-1,Ny+1,2:end-1],2), vmin=vmin, vmax=vmax, opacity=0.8)

  # surfaces 5,6 (z=z[1], z=z[end]
  mayavis.mesh(X[:,:,floor(Nz/2.0)],Y[:,:,floor(Nz/2.0)],Z[:,:,floor(Nz/2.0)], 
    scalars=phi.value[2:end-1,2:end-1,floor(Nz/2.0)+1], vmin=vmin, vmax=vmax, opacity=0.8)
  mayavis.mesh(X[:,:,Nz],Y[:,:,Nz],Z[:,:,Nz], 
    scalars=phi.value[2:end-1,2:end-1,Nz+1], vmin=vmin, vmax=vmax, opacity=0.8)

  mayavis.colorbar()
  mayavis.axes()
  return mayavis.screenshot()

end
end  