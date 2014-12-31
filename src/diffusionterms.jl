# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 27 December, 2014
# ===============================


# ========================= DIFFUSION TERM ================================
function diffusionTerm(D::FaceValue)
d = D.domain.dimension
if d==1
  M = diffusionTerm1D(D)
elseif d==1.5
  M = diffusionTermCylindrical1D(D)
elseif d==2
  M, Mx, My = diffusionTerm2D(D)
elseif d==2.5
  M, Mx, My = diffusionTermCylindrical2D(D)
elseif d==2.8
  M, Mx, My = diffusionTermRadial2D(D)
elseif d==3
  M, Mx, My, Mz = diffusionTerm3D(D)
elseif d==3.2
  M, Mx, My, Mz = diffusionTermCylindrical3D(D)
end
M
end

# ======================== 1D CARTESIAN DIFFUSION =========================
function diffusionTerm1D(D::FaceValue)
# D is a face variable

# extract data from the mesh structure
G = D.domain.numbering
Nx = D.domain.numberofcells[1]
dx = D.domain.cellsize[1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the diffusion coefficient data
Dx = D.xvalue

# reassign the east, west for code readability
De = Dx[2:Nx+1]
Dw = Dx[1:Nx]

# calculate the coefficients for the internal cells
AE = reshape(De/(dx*dx),Nx)
AW = reshape(Dw/(dx*dx),Nx)
APx = reshape(-(De+Dw)/(dx*dx),Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)
end

# ======================== 2D CARTESIAN DIFFUSION =========================
function diffusionTerm2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
G = D.domain.numbering
Nx = D.domain.numberofcells[1]
Ny = D.domain.numberofcells[2]
dx = D.domain.cellsize[1]
dy = D.domain.cellsize[2]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# extract the diffusion coefficient data
Dx = D.xvalue
Dy = D.yvalue

# reassign the east, west for code readability
De = Dx[2:Nx+1,:]
Dw = Dx[1:Nx,:]
Dn = Dy[:,2:Ny+1]
Ds = Dy[:,1:Ny]

# calculate the coefficients for the internal cells
AE = reshape(De/(dx*dx),mnx)
AW = reshape(Dw/(dx*dx),mnx)
AN = reshape(Dn/(dy*dy),mny)
AS = reshape(Ds/(dy*dy),mny)
APx = reshape(-(De+Dw)/(dx*dx),mnx)
APy = reshape(-(Dn+Ds)/(dy*dy),mny)

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

# ======================== 3D CARTESIAN DIFFUSION =========================
function diffusionTerm3D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
G = D.domain.numbering
Nx = D.domain.numberofcells[1]
Ny = D.domain.numberofcells[2]
Nz = D.domain.numberofcells[3]
dx = D.domain.cellsize[1]
dy = D.domain.cellsize[2]
dz = D.domain.cellsize[3]

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
Dx = D.xvalue
Dy = D.yvalue
Dz = D.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
De = Dx[2:Nx+1,:,:]
Dw = Dx[1:Nx,:,:]
Dn = Dy[:,2:Ny+1,:]
Ds = Dy[:,1:Ny,:]
Df = Dz[:,:,2:Nz+1]
Db = Dz[:,:,1:Nz]

# calculate the coefficients for the internal cells
AE = reshape(De/(dx*dx),mnx)
AW = reshape(Dw/(dx*dx),mnx)
AN = reshape(Dn/(dy*dy),mny)
AS = reshape(Ds/(dy*dy),mny)
AF = reshape(Df/(dz*dz),mnz)
AB = reshape(Db/(dz*dz),mnz)
APx = reshape(-(De+Dw)/(dx*dx),mnx)
APy = reshape(-(Dn+Ds)/(dy*dy),mny)
APz = reshape(-(Df+Db)/(dz*dz),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3);
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); 
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); 
		reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); 
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); 
		reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); 
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); 
		reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)];
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



# ======================== 1D Cylindrical DIFFUSION =========================
function diffusionTermCylindrical1D(D::FaceValue)
# D is a face variable

# extract data from the mesh structure
G = D.domain.numbering
Nx = D.domain.numberofcells[1]
dx = D.domain.cellsize[1]
rp = D.domain.cellcenters.x
rf = D.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# extract the diffusion coefficient data
Dx = D.xvalue

# reassign the east, west for code readability
De = Dx[2:Nx+1]
Dw = Dx[1:Nx]
re = rf[2:Nx+1]
rw = rf[1:Nx]

# calculate the coefficients for the internal cells
AE = reshape(re.*De./(dx*dx*rp),Nx)
AW = reshape(rw.*Dw./(dx*dx*rp),Nx)
APx = reshape(-(re.*De+rw.*Dw)./(dx*dx*rp),Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] = repmat(rowx_index,3)
jjx[1:3*Nx] = [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] = [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)
end



# ======================== 2D RADIAL DIFFUSION =========================
function diffusionTermRadial2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
G = D.domain.numbering
Nr = D.domain.numberofcells[1]
Ntheta = D.domain.numberofcells[2]
dr = D.domain.cellsize[1]
dtheta = D.domain.cellsize[2]
rp = repmat(D.domain.cellcenters.x, 1, Ntheta)
rf = repmat(D.domain.facecenters.x, 1, Ntheta)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# extract the diffusion coefficient data
Dx = D.xvalue
Dy = D.yvalue

# reassign the east, west for code readability
De = Dx[2:Nr+1,:]
Dw = Dx[1:Nr,:]
Dn = Dy[:,2:Ntheta+1]
Ds = Dy[:,1:Ntheta]
re = rf[2:Nr+1,:]         
rw = rf[1:Nr,:]

# calculate the coefficients for the internal cells
AE = reshape(re.*De./(dr*dr*rp),mnx)
AW = reshape(rw.*Dw./(dr*dr*rp),mnx)
AN = reshape(Dn./(dtheta*dtheta*rp.*rp),mny)
AS = reshape(Ds./(dtheta*dtheta*rp.*rp),mny)
APx = reshape(-(re.*De+rw.*Dw)./(dr*dr*rp),mnx)
APy = reshape(-(Dn+Ds)./(dtheta*dtheta*rp.*rp),mny)

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



# ======================== 2D CYLINDRICAL DIFFUSION =========================
function diffusionTermCylindrical2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
G = D.domain.numbering
Nr = D.domain.numberofcells[1]
Nz = D.domain.numberofcells[2]
dr = D.domain.cellsize[1]
dz = D.domain.cellsize[2]
rp = repmat(D.domain.cellcenters.x, 1, Nz)
rf = repmat(D.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# extract the diffusion coefficient data
Dx = D.xvalue
Dy = D.yvalue

# reassign the east, west for code readability
De = Dx[2:Nr+1,:]
Dw = Dx[1:Nr,:]
Dn = Dy[:,2:Nz+1]
Ds = Dy[:,1:Nz]
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# calculate the coefficients for the internal cells
AE = reshape(re.*De./(dr*dr*rp),mnx)
AW = reshape(rw.*Dw./(dr*dr*rp),mnx)
AN = reshape(Dn/(dz*dz),mny)
AS = reshape(Ds/(dz*dz),mny)
APx = reshape(-(re.*De+rw.*Dw)./(dr*dr*rp),mnx)
APy = reshape(-(Dn+Ds)/(dz*dz),mny)

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



# ======================== 3D CYLINDRICAL DIFFUSION =========================
function diffusionTermCylindrical3D(D::FaceValue)
# extract data from the mesh structure
G = D.domain.numbering
Nr = D.domain.numberofcells[1]
Ntheta = D.domain.numberofcells[2]
Nz = D.domain.numberofcells[3]
dr = D.domain.cellsize[1]
dtheta = D.domain.cellsize[2]
dz = D.domain.cellsize[3]
#rp = repmat(D.domain.cellcenters.x, Ntheta, Nz)
#rf = repmat(D.domain.facecenters.x, Ntheta, Nz)
rp = D.domain.cellcenters.x # use broadcasting
rf = D.domain.facecenters.x

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
Dx = D.xvalue
Dy = D.yvalue
Dz = D.zvalue

# reassign the east, west, north, and south velocity vectors for the 
# code readability
De = Dx[2:Nr+1,:,:]
Dw = Dx[1:Nr,:,:]
Dn = Dy[:,2:Ntheta+1,:]
Ds = Dy[:,1:Ntheta,:]
Df = Dz[:,:,2:Nz+1]
Db = Dz[:,:,1:Nz]
#re = rf[2:Nr+1,:,:]
#rw = rf[1:Nr,:,:]
re = rf[2:Nr+1] # use broadcasting
rw = rf[1:Nr]

# calculate the coefficients for the internal cells
AE = reshape(re.*De./(dr*dr*rp),mnx)
AW = reshape(rw.*Dw./(dr*dr*rp),mnx)
AN = reshape(Dn./(dtheta*dtheta*rp.*rp),mny)
AS = reshape(Ds./(dtheta*dtheta*rp.*rp),mny)
AF = reshape(Df/(dz*dz),mnz)
AB = reshape(Db/(dz*dz),mnz)
APx = reshape(-(re.*De+rw.*Dw)./(dr*dr*rp),mnx)
APy = reshape(-(Dn+Ds)./(dtheta*dtheta*rp.*rp),mny)
APz = reshape(-(Df+Db)/(dz*dz),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] = repmat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] = repmat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] = repmat(rowz_index,3)
jjx[1:3*mnx] = [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); 
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); 
		reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] = [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); 
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); 
		reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] = [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); 
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); 
		reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
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