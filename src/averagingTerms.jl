# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 27 December, 2014
# ===============================

# =========================================================
# Last changes:
# 27 December, 2014
#     - Radial and 3D cylindrical support
#     - improved averaging for cells witout a ghost value
#     - finished up the tvd averaging scheme
# =========================================================

# ================== Aritmetic averaging scheme ==================
function arithmeticMean(phi::CellValue)
# calculates the average values of a cell variable. The output is a
# face variable
d=phi.domain.dimension
if d==1 || d==1.5
  FaceValue(phi.domain,
    0.5*(phi.value[1:end-1]+phi.value[2:end]),
    [1.0],
    [1.0])
elseif d==2 || d==2.5 || d==2.8
  FaceValue(phi.domain,
    0.5*(phi.value[1:end-1,2:end-1]+phi.value[2:end,2:end-1]),
    0.5*(phi.value[2:end-1,1:end-1]+phi.value[2:end-1,2:end]),
    [1.0])
elseif d==3 || d==3.2
  FaceValue(phi.domain,
    0.5*(phi.value[1:end-1,2:end-1,2:end-1]+phi.value[2:end,2:end-1,2:end-1]),
    0.5*(phi.value[2:end-1,1:end-1,2:end-1]+phi.value[2:end-1,2:end,2:end-1]),
    0.5*(phi.value[2:end-1,2:end-1,1:end-1]+phi.value[2:end-1,2:end-1,2:end]))
end
end

# ================== Geometric averaging scheme ==================
function geometricMean(phi::CellValue)
# calculates the average values of a cell variable. The output is a
# face variable
d=phi.domain.dimension
if d==1 || d==1.5
  FaceValue(phi.domain,
    sqrt(phi.value[1:end-1].*phi.value[2:end]),
    [1.0],
    [1.0])
elseif d==2 || d==2.5 || d==2.8
  FaceValue(phi.domain,
    sqrt(phi.value[1:end-1,2:end-1].*phi.value[2:end,2:end-1]),
    sqrt(phi.value[2:end-1,1:end-1].*phi.value[2:end-1,2:end]),
    [1.0])
elseif d==3 || d==3.2
  FaceValue(phi.domain,
    sqrt(phi.value[1:end-1,2:end-1,2:end-1].*phi.value[2:end,2:end-1,2:end-1]),
    sqrt(phi.value[2:end-1,1:end-1,2:end-1].*phi.value[2:end-1,2:end,2:end-1]),
    sqrt(phi.value[2:end-1,2:end-1,1:end-1].*phi.value[2:end-1,2:end-1,2:end]))
end
end

# ================== Harmonic averaging scheme ==================
function harmonicMean(phi::CellValue)
# calculates the average values of a cell variable. The output is a
# face variable
d=phi.domain.dimension
if d==1 || d==1.5
  FaceValue(phi.domain,
    2.0./(1.0./phi.value[1:end-1]+1.0./phi.value[2:end]),
    [1.0],
    [1.0])
elseif d==2 || d==2.5 || d==2.8
  FaceValue(phi.domain,
    2.0./(1.0./phi.value[1:end-1,2:end-1]+1.0./phi.value[2:end,2:end-1]),
    2.0./(1.0./phi.value[2:end-1,1:end-1]+1.0./phi.value[2:end-1,2:end]),
    [1.0])
elseif d==3 || d==3.2
  FaceValue(phi.domain,
    2.0./(1.0./phi.value[1:end-1,2:end-1,2:end-1]+1.0./phi.value[2:end,2:end-1,2:end-1]),
    2.0./(1.0./phi.value[2:end-1,1:end-1,2:end-1]+1.0./phi.value[2:end-1,2:end,2:end-1]),
    2.0./(1.0./phi.value[2:end-1,2:end-1,1:end-1]+1.0./phi.value[2:end-1,2:end-1,2:end]))
end
end

# ================== Upwind averaging scheme ==================
function upwindMean(phi::CellValue, u::FaceValue)
# calculates the average values of a cell variable. The output is a
# face variable
d=phi.domain.dimension
phi_tmp = phi.value
if d==1 || d==1.5
  ux=u.xvalue
  # assign the value of the left boundary to the left ghost cell
  phi_tmp[1] = 0.5*(phi[1]+phi[2])
  # assign the value of the right boundary to the right ghost cell
  phi_tmp[end] = 0.5*(phi[end]+phi[end-1])
  FaceValue(phi.domain,
    (ux.>0.0).*phi_tmp[1:end-1]+(ux.<0.0).*phi_tmp[2:end]+
    0.5*(ux.==0).*(phi.value[1:end-1]+phi.value[2:end]),
    [1.0],
    [1.0])
elseif d==2 || d==2.5 || d==2.8
  ux = u.xvalue
  uy = u.yvalue
  # assign the value of the left boundary to the left ghost cells
  phi_tmp[1,:] = 0.5*(phi[1,:]+phi[2,:])
  # assign the value of the right boundary to the right ghost cells
  phi_tmp[end,:] = 0.5*(phi[end,:]+phi[end-1,:])
  # assign the value of the bottom boundary to the bottom ghost cells
  phi_tmp[:,1] = 0.5*(phi[:,1]+phi[:,2])
  # assign the value of the top boundary to the top ghost cells
  phi_tmp[:,end] = 0.5*(phi[:,end]+phi[:,end-1])
  FaceValue(phi.domain,
    (ux.>0.0).*phi_tmp[1:end-1,2:end-1]+
    (ux.<0.0).*phi_tmp[2:end,2:end-1]+
    0.5*(ux.==0.0).*(phi.value[1:end-1,2:end-1]+phi.value[2:end,2:end-1]),
    (uy.>0.0).*phi_tmp[2:end-1,1:end-1]+
    (uy.<0.0).*phi_tmp[2:end-1,2:end]+
    0.5*(uy.==0.0).*(phi.value[2:end-1,1:end-1]+phi.value[2:end-1,2:end]),
    [1.0])
elseif d==3 || d==3.2
  ux = u.xvalue
  uy = u.yvalue
  uz = u.zvalue
  # assign the value of the left boundary to the left ghost cells
  phi_tmp[1,:,:] = 0.5*(phi[1,:,:]+phi[2,:,:])
  # assign the value of the right boundary to the right ghost cells
  phi_tmp[end,:,:] = 0.5*(phi[end,:,:]+phi[end-1,:,:])
  # assign the value of the bottom boundary to the bottom ghost cells
  phi_tmp[:,1,:] = 0.5*(phi[:,1,:]+phi[:,2,:])
  # assign the value of the top boundary to the top ghost cells
  phi_tmp[:,end,:] = 0.5*(phi[:,end,:]+phi[:,end-1,:])
  # assign the value of the back boundary to the back ghost cells
  phi_tmp[:,:,1] = 0.5*(phi[:,:,1]+phi[:,:,2])
  # assign the value of the front boundary to the front ghost cells
  phi_tmp[:,:,end] = 0.5*(phi[:,:,end]+phi[:,:,end-1])
  FaceValue(phi.domain,
    (ux.>0.0).*phi_tmp[1:end-1,2:end-1,2:end-1]+
    (ux.<0.0).*phi_tmp[2:end,2:end-1,2:end-1]+
    0.5*(ux.==0.0).*(phi.value[1:end-1,2:end-1,2:end-1]+phi.value[2:end,2:end-1,2:end-1]),
    (uy.>0.0).*phi_tmp[2:end-1,1:end-1,2:end-1]+
    (uy.<0.0).*phi_tmp[2:end-1,2:end,2:end-1]+
    0.5*(uy.==0.0).*(phi.value[2:end-1,1:end-1,2:end-1]+phi.value[2:end-1,2:end,2:end-1]),
    (uz.>0.0).*phi_tmp[2:end-1,2:end-1,1:end-1]+
    (uz.<0.0).*phi_tmp[2:end-1,2:end-1,2:end]+
    0.5*(uz.==0.0).*(phi.value[2:end-1,2:end-1,1:end-1]+phi.value[2:end-1,2:end-1,2:end]))
end
end

# ================== TVD averaging scheme ==================
function tvdMean(phi::CellValue, u::FaceValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

d=phi.domain.dimension

if d == 1 || d == 1.5
  # extract data from the mesh structure
  Nx = u.domain.numberofcells[1]
  phi_p = zeros(Float64, Nx+1)
  phi_m = zeros(Float64, Nx+1)

  # extract the velocity data
  ux = u.xvalue

  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  dphi_p = phi[2:Nx+2]-phi[1:Nx+1]
  rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
  phi_p[2:Nx+1] = phi[2:Nx+1]+0.5*FL(rp).*(phi[3:Nx+2]-phi[2:Nx+1])
  phi_p[1] = (phi[1]+phi[2])/2.0 # left boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
  phi_m[1:Nx] = phi[2:Nx+1]+0.5*FL(rm).*(phi[1:Nx]-phi[2:Nx+1])
  phi_m[Nx+1] = (phi[end]+phi[end-1])/2.0 # right boundary

  FaceValue(phi.domain,
      (ux.>0.0).*phi_p+(ux.<0.0).*phi_m+
      0.5*(ux.==0).*(phi.value[1:end-1]+phi.value[2:end]),
      [1.0],
      [1.0])
elseif d==2 || d==2.5 || d==2.8
  # extract data from the mesh structure
  Nx = u.domain.numberofcells[1]
  Ny = u.domain.numberofcells[2]
  phi_p = zeros(Float64, Nx+1)
  phi_m = zeros(Float64, Nx+1)
  phiX_p = zeros(Float64, Nx+1, Ny)
  phiX_m = zeros(Float64, Nx+1,Ny)
  phiY_p = zeros(Float64, Nx,Ny+1)
  phiY_m = zeros(Float64, Nx,Ny+1)
  
  # extract the velocity data
  ux = u.xvalue
  uy = u.yvalue
  
  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  # x direction
  dphiX_p = phi[2:Nx+2, 2:Ny+1]-phi[1:Nx+1, 2:Ny+1]
  rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
  phiX_p[2:Nx+1,:] = phi[2:Nx+1, 2:Ny+1]+0.5*FL(rX_p).*
  (phi[3:Nx+2,2:Ny+1]-phi[2:Nx+1, 2:Ny+1])
  phiX_p[1, :] = (phi[1, 2:Ny+1]+phi[2, 2:Ny+1])/2.0  # left boundary
  # y direction
  dphiY_p = phi[2:Nx+1, 2:Ny+2]-phi[2:Nx+1, 1:Ny+1]
  rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
  phiY_p[:,2:Ny+1] = phi[2:Nx+1, 2:Ny+1]+0.5*FL(rY_p).*
      (phi[2:Nx+1,3:Ny+2]-phi[2:Nx+1, 2:Ny+1])
  phiY_p[:,1] = (phi[2:Nx+1,1]+phi[2:Nx+1,2])/2.0  # Bottom boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  # x direction
  rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
  phiX_m[1:Nx,:] = phi[2:Nx+1, 2:Ny+1]+0.5*FL(rX_m).*
      (phi[1:Nx, 2:Ny+1]-phi[2:Nx+1, 2:Ny+1])
  phiX_m[Nx+1,:] = (phi[end, 2:Ny+1]+phi[end-1, 2:Ny+1])/2.0  # right boundary
  # y direction
  rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
  phiY_m[:,1:Ny] = phi[2:Nx+1, 2:Ny+1]+0.5*FL(rY_m).*
      (phi[2:Nx+1, 1:Ny]-phi[2:Nx+1, 2:Ny+1])
  phiY_m[:, Ny+1] = (phi[2:Nx+1, end]+phi[2:Nx+1, end-1])/2.0  # top boundary
  
  
  FaceValue(phi.domain,
      (ux.>0.0).*phi_p+(ux.<0.0).*phi_m+
        0.5*(ux.==0.0).*(phi.value[1:end-1]+phi.value[2:end]),
      [1.0],
      [1.0])
      
  FaceValue(phi.domain,
      (ux.>0.0).*phiX_p+(ux.<0.0).*phiX_m+
         0.5*(ux.==0.0).*(phi[1:Nx+1,2:Ny+1]+phi[2:Nx+2,2:Ny+1]),
      (uy.>0.0).*phiY_p+(uy.<0.0).*phiY_m+
         0.5*(uy.==0.0).*(phi[2:Nx+1,1:Ny+1]+phi[2:Nx+1,2:Ny+2]),
      [1.0])
      
elseif d==3 || d==3.2
  # extract data from the mesh structure
  Nx = u.domain.numberofcells[1]
  Ny = u.domain.numberofcells[2]
  Nz = u.domain.numberofcells[3]
  
  # extract the velocity data
  ux = u.xvalue
  uy = u.yvalue
  uz = u.zvalue
  
  # define the tvd face vectors
  phiX_p = zeros(Float64, Nx+1,Ny,Nz)
  phiX_m = zeros(Float64, Nx+1,Ny,Nz)
  phiY_p = zeros(Float64, Nx,Ny+1,Nz)
  phiY_m = zeros(Float64, Nx,Ny+1,Nz)
  phiZ_p = zeros(Float64, Nx,Ny,Nz+1)
  phiZ_m = zeros(Float64, Nx,Ny,Nz+1)

  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  # x direction
  dphiX_p = phi[2:Nx+2, 2:Ny+1, 2:Nz+1]-phi[1:Nx+1, 2:Ny+1, 2:Nz+1]
  rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
  phiX_p[2:Nx+1,:,:] = phi[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rX_p).*
      (phi[3:Nx+2,2:Ny+1,2:Nz+1]-phi[2:Nx+1,2:Ny+1,2:Nz+1])
  phiX_p[1,:,:] = (phi[1,2:Ny+1,2:Nz+1]+phi[2,2:Ny+1,2:Nz+1])/2.0  # left boundary
  # y direction
  dphiY_p = phi[2:Nx+1, 2:Ny+2, 2:Nz+1]-phi[2:Nx+1, 1:Ny+1, 2:Nz+1]
  rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
  phiY_p[:,2:Ny+1,:] = phi[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rY_p).*
      (phi[2:Nx+1,3:Ny+2,2:Nz+1]-phi[2:Nx+1, 2:Ny+1,2:Nz+1])
  phiY_p[:,1,:] = (phi[2:Nx+1,1,2:Nz+1]+phi[2:Nx+1,2,2:Nz+1])/2.0  # Bottom boundary
  # z direction
  dphiZ_p = phi[2:Nx+1, 2:Ny+1, 2:Nz+2]-phi[2:Nx+1, 2:Ny+1, 1:Nz+1]
  rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
  phiZ_p[:,:,2:Nz+1] = phi[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rZ_p).*
      (phi[2:Nx+1,2:Ny+1,3:Nz+2]-phi[2:Nx+1,2:Ny+1,2:Nz+1])
  phiZ_p[:,:,1] = (phi[2:Nx+1,2:Ny+1,1]+phi[2:Nx+1,2:Ny+1,2])/2.0  # Back boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  # x direction
  rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
  phiX_m[1:Nx,:,:] = phi[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rX_m).*
      (phi[1:Nx, 2:Ny+1, 2:Nz+1]-phi[2:Nx+1, 2:Ny+1, 2:Nz+1])
  phiX_m[Nx+1,:,:] = (phi[end,2:Ny+1,2:Nz+1]+phi[end-1,2:Ny+1,2:Nz+1])/2.0  # right boundary
  # y direction
  rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
  phiY_m[:,1:Ny,:] = phi[2:Nx+1,2:Ny+1,2:Nz+1]+0.5*FL(rY_m).*
      (phi[2:Nx+1,1:Ny,2:Nz+1]-phi[2:Nx+1,2:Ny+1,2:Nz+1])
  phiY_m[:,Ny+1,:] = (phi[2:Nx+1, end,2:Nz+1]+phi[2:Nx+1, end-1,2:Nz+1])/2.0  # top boundary
  # z direction
  rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
  phiZ_m[:,:,1:Nz] = phi[2:Nx+1,2:Ny+1,2:Nz+1]+0.5*FL(rZ_m).*
      (phi[2:Nx+1,2:Ny+1,1:Nz]-phi[2:Nx+1,2:Ny+1,2:Nz+1])
  phiZ_m[:,:,Nz+1] = (phi[2:Nx+1,2:Ny+1,end]+phi[2:Nx+1,2:Ny+1,end-1])/2.0  # front boundary
  
  FaceValue(phi.domain,
      (ux.>0.0).*phiX_p+(ux.<0.0).*phiX_m+
         0.5*(ux.==0.0).*(phi[1:Nx+1,2:Ny+1,2:Nz+1]+phi[2:Nx+2,2:Ny+1,2:Nz+1]),
      (uy.>0.0).*phiY_p+(uy.<0).*phiY_m+
         0.5*(uy.==0.0).*(phi[2:Nx+1,1:Ny+1,2:Nz+1]+phi[2:Nx+1,2:Ny+2,2:Nz+1]),
      (uz.>0.0).*phiZ_p+(uz.<0).*phiZ_m+
         0.5*(uz.==0.0).*(phi[2:Nx+1,2:Ny+1,1:Nz+1]+phi[2:Nx+1,2:Ny+1,2:Nz+2]))
  
  
end

end