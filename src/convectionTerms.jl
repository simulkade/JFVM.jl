# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 29 December, 2014
# ===============================

# ======================================================
# 2014-12-29
#    - Import all the convection terms from matlab code
# ======================================================



function convectionTerm(u::FaceValue)
d = u.domain.dimension
if d==1
  M = convectionTerm1D(u)
elseif d==1.5
  M = convectionTermCylindrical1D(u)
elseif d==2
  M, Mx, My = convectionTerm2D(u)
elseif d==2.5
  M, Mx, My = convectionTermCylindrical2D(u)
elseif d==2.8
  M, Mx, My = convectionTermRadial2D(u)
elseif d==3
  M, Mx, My, Mz = convectionTerm3D(u)
elseif d==3.2
  M, Mx, My, Mz = convectionTermCylindrical3D(u)
end
M
end

function convectionUpwindTerm(u::FaceValue)
d = u.domain.dimension
if d==1
  M = convectionUpwindTerm1D(u)
elseif d==1.5
  M = convectionUpwindTermCylindrical1D(u)
elseif d==2
  M, Mx, My = convectionUpwindTerm2D(u)
elseif d==2.5
  M, Mx, My = convectionUpwindTermCylindrical2D(u)
elseif d==2.8
  M, Mx, My = convectionUpwindTermRadial2D(u)
elseif d==3
  M, Mx, My, Mz = convectionUpwindTerm3D(u)
elseif d==3.2
  M, Mx, My, Mz = convectionUpwindTermCylindrical3D(u)
end
M
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function)
d = u.domain.dimension
if d==1
  M, RHS = convectionTvdTerm1D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==1.5
  M, RHS = convectionTvdTermCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==2
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTerm2D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==2.5
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTermCylindrical2D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==2.8
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTermRadial2D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==3
  M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz = convectionTvdTerm3D(u::FaceValue, phi::CellValue, FL::Function)
elseif d==3.2
  M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz = convectionTvdTermCylindrical3D(u::FaceValue, phi::CellValue, FL::Function)
end
(M, RHS)
end

# =================== 1D Convection Terms =======================
# =============== 1D Convection Terms Cartesian Central =================
function convectionTerm1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]

# calculate the coefficients for the internal cells
AE = reshape(ue/(2.0*dx),Nx)
AW = reshape(-uw/(2.0*dx),Nx)
APx = reshape((ue-uw)/(2.0*dx),Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end


# =============== 1D Convection Terms Cylindrical Central =================
function convectionTermCylindrical1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]
re = rf[2:Nx+1]     
rw = rf[1:Nx]

# calculate the coefficients for the internal cells
AE = reshape(re.*ue./(2.0*dx*rp),Nx)
AW = reshape(-rw.*uw./(2.0*dx*rp),Nx)
APx = reshape((re.*ue-rw.*uw)./(2.0*dx.*rp),Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end


# =============== 1D Convection Term Cartesian Upwind =================
function convectionUpwindTerm1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)

# calculate the coefficients for the internal cells
AE = reshape(ue_min/dx,Nx)
AW = reshape(-uw_max/dx,Nx)
APx = reshape((ue_max-uw_min)/dx,Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1]-uw_max[1]/(2.0*dx)   
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0 
APx[end] = APx[end] + ue_min[end]/(2.0*dx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end

# =============== 1D Convection Terms Cylindrical Upwind =================
function convectionUpwindTermCylindrical1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]
re = rf[2:Nx+1]     
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)

# calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(dx*rp),Nx)
AW = reshape(-rw.*uw_max./(dx*rp),Nx)
APx = reshape((re.*ue_max-rw.*uw_min)./(dx*rp),Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1]-rw[1]*uw_max[1]/(2.0*dx*rp[1])   
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0 
APx[end] = APx[end] + re[end]*ue_min[end]/(2.0*dx*rp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end


# ================= 1D Convection Terms Cylindrical TVD =====================
function convectionTvdTermCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)
r = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue
# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = phi.value[2:Nx+2]-phi.value[1:Nx+1]
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] = 0.5*FL(rp).*(phi.value[3:Nx+2]-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] = 0.5*FL(rm).*(phi.value[1:Nx]-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]
re = rf[2:Nx+1]     
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)

# calculate the TVD correction term
RHS[2:Nx+1] = -(1.0./(dx*r)).*(re.*(ue_max.*psi_p[2:Nx+1]+ue_min.*psi_m[2:Nx+1])-
              rw.*(uw_max.*psi_p[1:Nx]+uw_min.*psi_m[1:Nx]))
              
# calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(dx*r),Nx)
AW = reshape(-rw.*uw_max./(dx*r),Nx)
APx = reshape((re.*ue_max-rw.*uw_min)./(dx*r),Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1]-rw[1]*uw_max[1]/(2.0*dx*r[1])   
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0 
APx[end] = APx[end] + re[end]*ue_min[end]/(2.0*dx*r[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

(M, RHS)
end


# ================= 1D Convection Terms Cartesian TVD =====================
function convectionTvdTerm1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
dx = u.domain.cellsize[1]
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the velocity data
ux = u.xvalue
# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = phi.value[2:Nx+2]-phi.value[1:Nx+1]
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] = 0.5*FL(rp).*(phi.value[3:Nx+2]-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] = 0.5*FL(rm).*(phi.value[1:Nx]-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# reassign the east, west for code readability
ue = ux[2:Nx+1]
uw = ux[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)

# calculate the TVD correction term
RHS[2:Nx+1] = -(1.0/dx)*((ue_max.*psi_p[2:Nx+1]+ue_min.*psi_m[2:Nx+1])-
              (uw_max.*psi_p[1:Nx]+uw_min.*psi_m[1:Nx]))
              
# calculate the coefficients for the internal cells
AE = reshape(ue_min/dx,Nx)
AW = reshape(-uw_max/dx,Nx)
APx = reshape((ue_max-uw_min)/dx,Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1]-uw_max[1]/(2.0*dx)   
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0 
APx[end] = APx[end] + ue_min[end]/(2.0*dx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

(M, RHS)
end



# =================== 2D Convection Terms =======================
# ================ 2D Convection Term Central ===================
function convectionTerm2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1,:]
uw = ux[1:Nx,:]
vn = uy[:,2:Ny+1]
vs = uy[:,1:Ny]

# calculate the coefficients for the internal cells
AE = reshape(ue/(2.0*dx),mnx)
AW = reshape(-uw/(2.0*dx),mnx)
AN = reshape(vn/(2.0*dy),mny)
AS = reshape(-vs/(2.0*dy),mny)
APx = reshape((ue-uw)/(2.0*dx),mnx)
APy = reshape((vn-vs)/(2.0*dy),mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term Upwind ===================
function convectionUpwindTerm2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nx+1,:]
uw = ux[1:Nx,:]
vn = uy[:,2:Ny+1]
vs = uy[:,1:Ny]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, not reshape
AE = ue_min/dx
AW = -uw_max/dx
AN = vn_min/dy
AS = -vs_max/dy
APx = (ue_max-uw_min)/dx
APy = (vn_max-vs_min)/dy

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-uw_max[1,:]/(2.0*dx)
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+ue_min[end,:]/(2.0*dx)
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]/(2.0*dy)
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]/(2.0*dy)

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term TVD =======================
function convectionTvdTerm2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]
psiX_p = zeros(Nx+1,Ny)
psiX_m = zeros(Nx+1,Ny)
psiY_p = zeros(Nx,Ny+1)
psiY_m = zeros(Nx,Ny+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = phi.value[2:Nx+2, 2:Ny+1]-phi.value[1:Nx+1, 2:Ny+1]
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nx+1,:] = 0.5*FL(rX_p).*(phi.value[3:Nx+2,2:Ny+1]-
		    phi.value[2:Nx+1, 2:Ny+1])
psiX_p[1, :] = 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = phi.value[2:Nx+1, 2:Ny+2]-phi.value[2:Nx+1, 1:Ny+1]
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ny+1] = 0.5*FL(rY_p).*(phi.value[2:Nx+1,3:Ny+2]-
		  phi.value[2:Nx+1, 2:Ny+1])
psiY_p[:,1] = 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nx,:] = 0.5*FL(rX_m).*(phi.value[1:Nx, 2:Ny+1]-
		phi.value[2:Nx+1, 2:Ny+1])
psiX_m[Nx+1,:] = 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ny] = 0.5*FL(rY_m).*(phi.value[2:Nx+1, 1:Ny]-
	      phi.value[2:Nx+1, 2:Ny+1])
psiY_m[:, Ny+1] = 0.0 # top boundary will be handled in the main matrix

# reassign the east, west for code readability
ue = ux[2:Nx+1,:]
uw = ux[1:Nx,:]
vn = uy[:,2:Ny+1]
vs = uy[:,1:Ny]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, not reshape
AE = ue_min/dx
AW = -uw_max/dx
AN = vn_min/dy
AS = -vs_max/dy
APx = (ue_max-uw_min)/dx
APy = (vn_max-vs_min)/dy

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-uw_max[1,:]/(2.0*dx)
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+ue_min[end,:]/(2.0*dx)
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]/(2.0*dy)
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]/(2.0*dy)

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1/dx)*((ue_max.*psiX_p[2:Nx+1,:]+ue_min.*psiX_m[2:Nx+1,:])-
              (uw_max.*psiX_p[1:Nx,:]+uw_min.*psiX_m[1:Nx,:]))
div_y = -(1/dy)*((vn_max.*psiY_p[:,2:Ny+1]+vn_min.*psiY_m[:,2:Ny+1])-
              (vs_max.*psiY_p[:,1:Ny]+vs_min.*psiY_m[:,1:Ny]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2))
RHSx = zeros((Nx+2)*(Ny+2))
RHSy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
RHS[rowx_index] = reshape(div_x+div_y,Nx*Ny)
RHSx[rowx_index] = reshape(div_x,Nx*Ny)
RHSy[rowy_index] = reshape(div_y,Nx*Ny)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end




# ================ 2D Convection Term Cylindrical Central ===================
function convectionTermCylindrical2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Nz = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dz = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Nz)
rf = repmat(u.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Nz+1]
vs = uy[:,1:Nz]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# calculate the coefficients for the internal cells
AE = reshape(re.*ue./(2.0*dr*rp),mnx)
AW = reshape(-rw.*uw./(2.0*dr*rp),mnx)
AN = reshape(vn/(2.0*dz),mny)
AS = reshape(-vs/(2.0*dz),mny)
APx = reshape((re.*ue-rw.*uw)./(2.0*dr*rp),mnx)
APy = reshape((vn-vs)/(2.0*dz),mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, Mx, My)
end



# ================ 2D Convection Term Cylindrical Upwind ==============================
function convectionUpwindTermCylindrical2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Nz = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dz = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Nz)
rf = repmat(u.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Nz+1]
vs = uy[:,1:Nz]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min/dz
AS = -vs_max/dz
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-rw[1,:].*uw_max[1,:]./(2.0*dr*rp[1,:])
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+re[end,:].*ue_min[end,:]./(2.0*dr*rp[end,:])
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]/(2.0*dz)
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]/(2.0*dz)

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term Cylindrical TVD =======================
function convectionTvdTermCylindrical2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Nz = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dz = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Nz)
rf = repmat(u.domain.facecenters.x, 1, Nz)
psiX_p = zeros(Nr+1,Nz)
psiX_m = zeros(Nr+1,Nz)
psiY_p = zeros(Nr,Nz+1)
psiY_m = zeros(Nr,Nz+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = phi.value[2:Nr+2, 2:Nz+1]-phi.value[1:Nr+1, 2:Nz+1]
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] = 0.5*FL(rX_p).*(phi.value[3:Nr+2,2:Nz+1]-
		    phi.value[2:Nr+1, 2:Nz+1])
psiX_p[1, :] = 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = phi.value[2:Nr+1, 2:Nz+2]-phi.value[2:Nr+1, 1:Nz+1]
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Nz+1] = 0.5*FL(rY_p).*(phi.value[2:Nr+1,3:Nz+2]-
		  phi.value[2:Nr+1, 2:Nz+1])
psiY_p[:,1] = 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] = 0.5*FL(rX_m).*(phi.value[1:Nr, 2:Nz+1]-
		phi.value[2:Nr+1, 2:Nz+1])
psiX_m[Nr+1,:] = 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Nz] = 0.5*FL(rY_m).*(phi.value[2:Nr+1, 1:Nz]-
	      phi.value[2:Nr+1, 2:Nz+1])
psiY_m[:, Nz+1] = 0.0 # top boundary will be handled in the main matrix

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Nz+1]
vs = uy[:,1:Nz]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, not reshape
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min/dz
AS = -vs_max/dz
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-rw[1,:].*uw_max[1,:]./(2.0*dr*rp[1,:])
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+re[end,:].*ue_min[end,:]./(2.0*dr*rp[end,:])
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]/(2.0*dz)
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]/(2.0*dz)

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1./(dr*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:]+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:]+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1/dz)*((vn_max.*psiY_p[:,2:Nz+1]+vn_min.*psiY_m[:,2:Nz+1])-
              (vs_max.*psiY_p[:,1:Nz]+vs_min.*psiY_m[:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
RHS[rowx_index] = reshape(div_x+div_y,Nr*Nz)
RHSx[rowx_index] = reshape(div_x,Nr*Nz)
RHSy[rowy_index] = reshape(div_y,Nr*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end


# ================ 2D Convection Term Radial Central ===================
function convectionTermRadial2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Ntheta)
rf = repmat(u.domain.facecenters.x, 1, Ntheta)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Ntheta+1]
vs = uy[:,1:Ntheta]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# calculate the coefficients for the internal cells
AE = reshape(re.*ue./(2.0*dr*rp),mnx)
AW = reshape(-rw.*uw./(2.0*dr*rp),mnx)
AN = reshape(vn./(2.0*dtheta*rp),mny)
AS = reshape(-vs./(2.0*dtheta*rp),mny)
APx = reshape((re.*ue-rw.*uw)./(2.0*dr*rp),mnx)
APy = reshape((vn-vs)./(2.0*dtheta*rp),mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, Mx, My)
end



# ================ 2D Convection Term Radial Upwind ==============================
function convectionUpwindTermRadial2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Ntheta)
rf = repmat(u.domain.facecenters.x, 1, Ntheta)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Ntheta+1]
vs = uy[:,1:Ntheta]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min./(dtheta*rp)
AS = -vs_max./(dtheta*rp)
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)./(dtheta*rp)

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-rw[1,:].*uw_max[1,:]./(2.0*dr*rp[1,:])
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+re[end,:].*ue_min[end,:]./(2.0*dr*rp[end,:])
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]./(2.0*dtheta*rp[:,1])
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]./(2.0*dtheta*rp[:,end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term Radial TVD =======================
function convectionTvdTermRadial2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
rp = repmat(u.domain.cellcenters.x, 1, Ntheta)
rf = repmat(u.domain.facecenters.x, 1, Ntheta)
psiX_p = zeros(Nr+1,Ntheta)
psiX_m = zeros(Nr+1,Ntheta)
psiY_p = zeros(Nr,Ntheta+1)
psiY_m = zeros(Nr,Ntheta+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# extract the velocity data
ux = u.xvalue
uy = u.yvalue

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = phi.value[2:Nr+2, 2:Ntheta+1]-phi.value[1:Nr+1, 2:Ntheta+1]
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] = 0.5*FL(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1]-
		    phi.value[2:Nr+1, 2:Ntheta+1])
psiX_p[1, :] = 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = phi.value[2:Nr+1, 2:Ntheta+2]-phi.value[2:Nr+1, 1:Ntheta+1]
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ntheta+1] = 0.5*FL(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2]-
		  phi.value[2:Nr+1, 2:Ntheta+1])
psiY_p[:,1] = 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] = 0.5*FL(rX_m).*(phi.value[1:Nr, 2:Ntheta+1]-
		phi.value[2:Nr+1, 2:Ntheta+1])
psiX_m[Nr+1,:] = 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ntheta] = 0.5*FL(rY_m).*(phi.value[2:Nr+1, 1:Ntheta]-
	      phi.value[2:Nr+1, 2:Ntheta+1])
psiY_m[:, Ntheta+1] = 0.0 # top boundary will be handled in the main matrix

# reassign the east, west for code readability
ue = ux[2:Nr+1,:]
uw = ux[1:Nr,:]
vn = uy[:,2:Ntheta+1]
vs = uy[:,1:Ntheta]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0.0)
ue_max = max(ue,0.0)
uw_min = min(uw,0.0)
uw_max = max(uw,0.0)
vn_min = min(vn,0.0)
vn_max = max(vn,0.0)
vs_min = min(vs,0.0)
vs_max = max(vs,0.0)

# calculate the coefficients for the internal cells, not reshape
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min./(dtheta*rp)
AS = -vs_max./(dtheta*rp)
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)./(dtheta*rp)

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] = APx[1,:]-rw[1,:].*uw_max[1,:]./(2.0*dr*rp[1,:])
AW[1,:] = AW[1,:]/2.0
# Right boundary:
AE[end,:] = AE[end,:]/2.0
APx[end,:] = APx[end,:]+re[end,:].*ue_min[end,:]./(2.0*dr*rp[end,:])
# Bottom boundary:
APy[:,1] = APy[:,1]-vs_max[:,1]./(2.0*dtheta*rp[:,1])
AS[:,1] = AS[:,1]/2.0
# Top boundary:
AN[:,end] = AN[:,end]/2.0
APy[:,end] = APy[:,end]+vn_min[:,end]./(2.0*dtheta*rp[:,end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1./(dr*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:]+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:]+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1./(dtheta*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1]+vn_min.*psiY_m[:,2:Ntheta+1])-
              (vs_max.*psiY_p[:,1:Ntheta]+vs_min.*psiY_m[:,1:Ntheta]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2))
RHSx = zeros((Nr+2)*(Ntheta+2))
RHSy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
RHS[rowx_index] = reshape(div_x+div_y,Nr*Ntheta)
RHSx[rowx_index] = reshape(div_x,Nr*Ntheta)
RHSy[rowy_index] = reshape(div_y,Nr*Ntheta)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end


# ============================= 3D Convection Term ===============================
function convectionTerm3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]
dz = u.domain.cellsize[3]

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nx+1,:,:]
uw = ux[1:Nx,:,:]
vn = uy[:,2:Ny+1,:]
vs = uy[:,1:Ny,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]

# calculate the coefficients for the internal cells
AE = reshape(ue/(2*dx),mnx)
AW = reshape(-uw/(2*dx),mnx)
AN = reshape(vn/(2*dy),mny)
AS = reshape(-vs/(2*dy),mny)
AF = reshape(wf/(2*dz),mnz)
AB = reshape(-wb/(2*dz),mnz)
APx = reshape((ue-uw)/(2*dx),mnx)
APy = reshape((vn-vs)/(2*dy),mny)
APz = reshape((wf-wb)/(2*dz),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz

(M, Mx, My, Mz)
end


# ============================= 3D Convection Upwind Term ===============================
function convectionUpwindTerm3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]
dz = u.domain.cellsize[3]

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nx+1,:,:]
uw = ux[1:Nx,:,:]
vn = uy[:,2:Ny+1,:]
vs = uy[:,1:Ny,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0)
ue_max = max(ue,0)
uw_min = min(uw,0)
uw_max = max(uw,0)
vn_min = min(vn,0)
vn_max = max(vn,0)
vs_min = min(vs,0)
vs_max = max(vs,0)
wf_min = min(wf,0)
wf_max = max(wf,0)
wb_min = min(wb,0)
wb_max = max(wb,0)

# calculate the coefficients for the internal cells
AE = ue_min/dx
AW = -uw_max/dx
AN = vn_min/dy
AS = -vs_max/dy
AF = wf_min/dz
AB = -wb_max/dz
APx = (ue_max-uw_min)/dx
APy = (vn_max-vs_min)/dy
APz = (wf_max-wb_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] = APx[1,:,:]-uw_max[1,:,:]/(2.0*dx)
AW[1,:,:] = AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] = AE[end,:,:]/2.0
APx[end,:,:] = APx[end,:,:]+ue_min[end,:,:]/(2.0*dx)
# Bottom boundary:
APy[:,1,:] = APy[:,1,:]-vs_max[:,1,:]/(2.0*dy)
AS[:,1,:] = AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] = AN[:,end,:]/2.0
APy[:,end,:] = APy[:,end,:]+vn_min[:,end,:]/(2.0*dy)
# Back boundary:
APz[:,:,1] = APz[:,:,1]-wb_max[:,:,1]/(2.0*dz)
AB[:,:,1] = AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] = AF[:,:,end]/2.0
APz[:,:,end] = APz[:,:,end]+wf_min[:,:,end]/(2.0*dz)

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end



# ============================= 3D Convection TVD Term ===============================
function convectionTvdTerm3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)
# extract data from the mesh structure
G = u.domain.numbering
Nx = u.domain.numberofcells[1]
Ny = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dx = u.domain.cellsize[1]
dy = u.domain.cellsize[2]
dz = u.domain.cellsize[3]
psiX_p = zeros(Nx+1,Ny,Nz)
psiX_m = zeros(Nx+1,Ny,Nz)
psiY_p = zeros(Nx,Ny+1,Nz)
psiY_m = zeros(Nx,Ny+1,Nz)
psiZ_p = zeros(Nx,Ny,Nz+1)
psiZ_m = zeros(Nx,Ny,Nz+1)

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = phi.value[2:Nx+2, 2:Ny+1, 2:Nz+1]-phi.value[1:Nx+1, 2:Ny+1, 2:Nz+1]
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nx+1,:,:] = 0.5*FL(rX_p).*(phi.value[3:Nx+2,2:Ny+1,2:Nz+1]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiX_p[1,:,:] = 0  # left boundary
# y direction
dphiY_p = phi.value[2:Nx+1, 2:Ny+2, 2:Nz+1]-phi.value[2:Nx+1, 1:Ny+1, 2:Nz+1]
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ny+1,:] = 0.5*FL(rY_p).*(phi.value[2:Nx+1,3:Ny+2,2:Nz+1]-phi.value[2:Nx+1, 2:Ny+1,2:Nz+1])
psiY_p[:,1,:] = 0.0  # Bottom boundary
# z direction
dphiZ_p = phi.value[2:Nx+1, 2:Ny+1, 2:Nz+2]-phi.value[2:Nx+1, 2:Ny+1, 1:Nz+1]
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] = 0.5*FL(rZ_p).*(phi.value[2:Nx+1,2:Ny+1,3:Nz+2]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_p[:,:,1] = 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nx,:,:] = 0.5*FL(rX_m).*(phi.value[1:Nx, 2:Ny+1, 2:Nz+1]-phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1])
psiX_m[Nx+1,:,:] = 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ny,:] = 0.5*FL(rY_m).*(phi.value[2:Nx+1,1:Ny,2:Nz+1]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiY_m[:,Ny+1,:] = 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] = 0.5*FL(rZ_m).*(phi.value[2:Nx+1,2:Ny+1,1:Nz]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_m[:,:,Nz+1] = 0.0  # front boundary

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nx+1,:,:]
uw = ux[1:Nx,:,:]
vn = uy[:,2:Ny+1,:]
vs = uy[:,1:Ny,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0)
ue_max = max(ue,0)
uw_min = min(uw,0)
uw_max = max(uw,0)
vn_min = min(vn,0)
vn_max = max(vn,0)
vs_min = min(vs,0)
vs_max = max(vs,0)
wf_min = min(wf,0)
wf_max = max(wf,0)
wb_min = min(wb,0)
wb_max = max(wb,0)

# calculate the coefficients for the internal cells
AE = ue_min/dx
AW = -uw_max/dx
AN = vn_min/dy
AS = -vs_max/dy
AF = wf_min/dz
AB = -wb_max/dz
APx = (ue_max-uw_min)/dx
APy = (vn_max-vs_min)/dy
APz = (wf_max-wb_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] = APx[1,:,:]-uw_max[1,:,:]/(2.0*dx)
AW[1,:,:] = AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] = AE[end,:,:]/2.0
APx[end,:,:] = APx[end,:,:]+ue_min[end,:,:]/(2.0*dx)
# Bottom boundary:
APy[:,1,:] = APy[:,1,:]-vs_max[:,1,:]/(2.0*dy)
AS[:,1,:] = AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] = AN[:,end,:]/2.0
APy[:,end,:] = APy[:,end,:]+vn_min[:,end,:]/(2.0*dy)
# Back boundary:
APz[:,:,1] = APz[:,:,1]-wb_max[:,:,1]/(2.0*dz)
AB[:,:,1] = AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] = AF[:,:,end]/2.0
APz[:,:,end] = APz[:,:,end]+wf_min[:,:,end]/(2.0*dz)

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# calculate the TVD correction term
div_x = -(1/dx)*((ue_max.*psiX_p[2:Nx+1,:,:]+ue_min.*psiX_m[2:Nx+1,:,:])-
              (uw_max.*psiX_p[1:Nx,:,:]+uw_min.*psiX_m[1:Nx,:,:]))
div_y = -(1/dy)*((vn_max.*psiY_p[:,2:Ny+1,:]+vn_min.*psiY_m[:,2:Ny+1,:])-
              (vs_max.*psiY_p[:,1:Ny,:]+vs_min.*psiY_m[:,1:Ny,:]))
div_z = -(1/dz)*((wf_max.*psiZ_p[:,:,2:Nz+1]+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz]+wb_min.*psiZ_m[:,:,1:Nz]))
          
# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
row_index = rowx_index
RHS[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSx[rowx_index] = reshape(div_x,Nx*Ny*Nz)
RHSy[rowy_index] = reshape(div_y,Nx*Ny*Nz)
RHSz[rowz_index] = reshape(div_z,Nx*Ny*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz;

(M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz)
end


# ============================= 3D Cylindrical Convection Term ===============================
function convectionTermCylindrical3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
dz = u.domain.cellsize[3]
#rp = repmat(u.domain.cellcenters.x, 1, Ntheta, Nz)
#rf = repmat(u.domain.facecenters.x, 1, Ntheta, Nz)
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nr+1,:,:]
uw = ux[1:Nr,:,:]
vn = uy[:,2:Ntheta+1,:]
vs = uy[:,1:Ntheta,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# calculate the coefficients for the internal cells
AE = reshape(re.*ue./(2.0*dr*rp),mnx)
AW = reshape(-rw.*uw./(2.0*dr*rp),mnx)
AN = reshape(vn./(2.0*dtheta*rp),mny)
AS = reshape(-vs./(2.0*dtheta*rp),mny)
AF = reshape(wf/(2.0*dz),mnz)
AB = reshape(-wb/(2.0*dz),mnz)
APx = reshape((re.*ue-rw.*uw)./(2.0*dr*rp),mnx)
APy = reshape((vn-vs)./(2.0*dtheta*rp),mny)
APz = reshape((wf-wb)/(2.0*dz),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz

(M, Mx, My, Mz)
end


# ============================= 3D Cylindrical Convection Upwind Term ===============================
function convectionUpwindTermCylindrical3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
dz = u.domain.cellsize[3]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nr+1,:,:]
uw = ux[1:Nr,:,:]
vn = uy[:,2:Ntheta+1,:]
vs = uy[:,1:Ntheta,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0)
ue_max = max(ue,0)
uw_min = min(uw,0)
uw_max = max(uw,0)
vn_min = min(vn,0)
vn_max = max(vn,0)
vs_min = min(vs,0)
vs_max = max(vs,0)
wf_min = min(wf,0)
wf_max = max(wf,0)
wb_min = min(wb,0)
wb_max = max(wb,0)

# calculate the coefficients for the internal cells
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min./(dtheta*rp)
AS = -vs_max./(dtheta*rp)
AF = wf_min/dz
AB = -wb_max/dz
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)./(dtheta*rp)
APz = (wf_max-wb_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] = APx[1,:,:]-rw[1,:,:].*uw_max[1,:,:]./(2.0*dr*rp[1,:,:])
AW[1,:,:] = AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] = AE[end,:,:]/2.0
APx[end,:,:] = APx[end,:,:]+re[end,:,:].*ue_min[end,:,:]./(2.0*dr*rp[end,:,:])
# Bottom boundary:
APy[:,1,:] = APy[:,1,:]-vs_max[:,1,:]./(2.0*dtheta*rp[:,1,:])
AS[:,1,:] = AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] = AN[:,end,:]/2.0
APy[:,end,:] = APy[:,end,:]+vn_min[:,end,:]./(2.0*dtheta*rp[:,end,:])
# Back boundary:
APz[:,:,1] = APz[:,:,1]-wb_max[:,:,1]/(2.0*dz)
AB[:,:,1] = AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] = AF[:,:,end]/2.0
APz[:,:,end] = APz[:,:,end]+wf_min[:,:,end]/(2.0*dz)

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end



# ============================= 3D Cylindrical Convection TVD Term ===============================
function convectionTvdTermCylindrical3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)
# extract data from the mesh structure
G = u.domain.numbering
Nr = u.domain.numberofcells[1]
Ntheta = u.domain.numberofcells[2]
Nz = u.domain.numberofcells[3]
dr = u.domain.cellsize[1]
dtheta = u.domain.cellsize[2]
dz = u.domain.cellsize[3]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x
psiX_p = zeros(Nr+1,Ntheta,Nz)
psiX_m = zeros(Nr+1,Ntheta,Nz)
psiY_p = zeros(Nr,Ntheta+1,Nz)
psiY_m = zeros(Nr,Ntheta+1,Nz)
psiZ_p = zeros(Nr,Ntheta,Nz+1)
psiZ_m = zeros(Nr,Ntheta,Nz+1)

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

# extract the velocity data 
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
ux = u.xvalue
uy = u.yvalue
uz = u.zvalue

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = phi.value[2:Nr+2, 2:Ntheta+1, 2:Nz+1]-phi.value[1:Nr+1, 2:Ntheta+1, 2:Nz+1]
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nr+1,:,:] = 0.5*FL(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1,2:Nz+1]-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiX_p[1,:,:] = 0  # left boundary
# y direction
dphiY_p = phi.value[2:Nr+1, 2:Ntheta+2, 2:Nz+1]-phi.value[2:Nr+1, 1:Ntheta+1, 2:Nz+1]
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ntheta+1,:] = 0.5*FL(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2,2:Nz+1]-phi.value[2:Nr+1, 2:Ntheta+1,2:Nz+1])
psiY_p[:,1,:] = 0.0  # Bottom boundary
# z direction
dphiZ_p = phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+2]-phi.value[2:Nr+1, 2:Ntheta+1, 1:Nz+1]
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] = 0.5*FL(rZ_p).*(phi.value[2:Nr+1,2:Ntheta+1,3:Nz+2]-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_p[:,:,1] = 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nr,:,:] = 0.5*FL(rX_m).*(phi.value[1:Nr, 2:Ntheta+1, 2:Nz+1]-phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+1])
psiX_m[Nr+1,:,:] = 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ntheta,:] = 0.5*FL(rY_m).*(phi.value[2:Nr+1,1:Ntheta,2:Nz+1]-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiY_m[:,Ntheta+1,:] = 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] = 0.5*FL(rZ_m).*(phi.value[2:Nr+1,2:Ntheta+1,1:Nz]-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_m[:,:,Nz+1] = 0.0  # front boundary

# reassign the east, west, north, and south velocity vectors for the 
# code readability
ue = ux[2:Nr+1,:,:]
uw = ux[1:Nr,:,:]
vn = uy[:,2:Ntheta+1,:]
vs = uy[:,1:Ntheta,:]
wf = uz[:,:,2:Nz+1]
wb = uz[:,:,2:Nz+1]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# find the velocity direction for the upwind scheme
ue_min = min(ue,0)
ue_max = max(ue,0)
uw_min = min(uw,0)
uw_max = max(uw,0)
vn_min = min(vn,0)
vn_max = max(vn,0)
vs_min = min(vs,0)
vs_max = max(vs,0)
wf_min = min(wf,0)
wf_max = max(wf,0)
wb_min = min(wb,0)
wb_max = max(wb,0)

# calculate the coefficients for the internal cells
AE = re.*ue_min./(dr*rp)
AW = -rw.*uw_max./(dr*rp)
AN = vn_min./(dtheta*rp)
AS = -vs_max./(dtheta*rp)
AF = wf_min/dz
AB = -wb_max/dz
APx = (re.*ue_max-rw.*uw_min)./(dr*rp)
APy = (vn_max-vs_min)./(dtheta*rp)
APz = (wf_max-wb_min)/dz

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] = APx[1,:,:]-rw[1,:,:].*uw_max[1,:,:]./(2.0*dr*rp[1,:,:])
AW[1,:,:] = AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] = AE[end,:,:]/2.0
APx[end,:,:] = APx[end,:,:]+re[end,:,:].*ue_min[end,:,:]./(2.0*dr*rp[end,:,:])
# Bottom boundary:
APy[:,1,:] = APy[:,1,:]-vs_max[:,1,:]./(2.0*dtheta*rp[:,1,:])
AS[:,1,:] = AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] = AN[:,end,:]/2.0
APy[:,end,:] = APy[:,end,:]+vn_min[:,end,:]./(2.0*dtheta*rp[:,end,:])
# Back boundary:
APz[:,:,1] = APz[:,:,1]-wb_max[:,:,1]/(2.0*dz)
AB[:,:,1] = AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] = AF[:,:,end]/2.0
APz[:,:,end] = APz[:,:,end]+wf_min[:,:,end]/(2.0*dz)

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] = [AW; APx; AE]
sy[1:3*mny] = [AS; APy; AN]
sz[1:3*mnz] = [AB; APz; AF]

# calculate the TVD correction term
div_x = -(1./(dr*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:,:]+ue_min.*psiX_m[2:Nr+1,:,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:,:]+uw_min.*psiX_m[1:Nr,:,:]))
div_y = -(1./(dtheta*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1,:]+vn_min.*psiY_m[:,2:Ntheta+1,:])-
              (vs_max.*psiY_p[:,1:Ntheta,:]+vs_min.*psiY_m[:,1:Ntheta,:]))
div_z = -(1/dz)*((wf_max.*psiZ_p[:,:,2:Nz+1]+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz]+wb_min.*psiZ_m[:,:,1:Nz]))
          
# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSz = zeros((Nr+2)*(Ntheta+2)*(Nz+2))

# assign the values of the RHS vector
row_index = rowx_index
RHS[row_index] = reshape(div_x+div_y+div_z,Nr*Ntheta*Nz)
RHSx[rowx_index] = reshape(div_x,Nr*Ntheta*Nz)
RHSy[rowy_index] = reshape(div_y,Nr*Ntheta*Nz)
RHSz[rowz_index] = reshape(div_z,Nr*Ntheta*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz;

(M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz)
end
