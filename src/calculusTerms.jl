# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 29 December, 2014
# ===============================

# ================================================================
# Changes:
#    2014-12-30 added 2D radial and 3D cylindrical grids
#    2014-12-31 many many bugs fixed
# ================================================================


function divergenceTerm(F::FaceValue)
d = F.domain.dimension
if d==1
  RHSdiv = divergenceTerm1D(F)
elseif d==1.5
  RHSdiv = divergenceTermCylindrical1D(F)
elseif d==2
  RHSdiv, RHSdivx, RHSdivy = divergenceTerm2D(F)
elseif d==2.5
  RHSdiv, RHSdivx, RHSdivy = divergenceTermCylindrical2D(F)
elseif d==2.8
  RHSdiv, RHSdivx, RHSdivy = divergenceTermRadial2D(F)
elseif d==3
  RHSdiv, RHSdivx, RHSdivy, RHSdivz = divergenceTerm3D(F)
elseif d==3.2
  RHSdiv, RHSdivx, RHSdivy, RHSdivz = divergenceTermCylindrical3D(F)
end
RHSdiv

end


# =============== Divergence 1D Term ============================
function divergenceTerm1D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nx = F.domain.numberofcells[1]
dx = F.domain.cellsize[1]
Fx = F.xvalue

# define the vector of cell index
row_index = reshape(G[2:Nx+1],Nx) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nx+1]
Fw = Fx[1:Nx]

# compute the divergence
div_x = (Fe-Fw)/dx

# define the RHS Vector
RHSdiv = zeros(Nx+2)

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x,Nx)

end


# =============== Divergence Cylindrical 1D Term ============================
function divergenceTermCylindrical1D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nx = F.domain.numberofcells[1]
dx = F.domain.cellsize[1]
Fx = F.xvalue
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x

# define the vector of cell index
row_index = reshape(G[2:Nx+1],Nx) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nx+1]
Fw = Fx[1:Nx]
re = rf[2:Nx+1]
rw = rf[1:Nx]

# compute the divergence
div_x = (re.*Fe-rw.*Fw)./(dx*rp)

# define the RHS Vector
RHSdiv = zeros(Nx+2)

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x,Nx)

end

# =============== Divergence 2D Term ============================
function divergenceTerm2D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nx = F.domain.numberofcells[1]
Ny = F.domain.numberofcells[2]
dx = F.domain.cellsize[1]
dy = F.domain.cellsize[2]
Fx = F.xvalue
Fy = F.yvalue

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nx+1,:]
Fw = Fx[1:Nx,:]
Fn = Fy[:,2:Ny+1]
Fs = Fy[:,1:Ny]

# compute the divergence
div_x = (Fe - Fw)/dx
div_y = (Fn - Fs)/dy

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2))
RHSdivx = zeros((Nx+2)*(Ny+2))
RHSdivy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nx*Ny)
RHSdivx[row_index] = reshape(div_x,Nx*Ny)
RHSdivy[row_index] = reshape(div_y,Nx*Ny)

(RHSdiv, RHSdivx, RHSdivy)

end



# =============== Divergence 2D Cylindrical Term ============================
function divergenceTermCylindrical2D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nr = F.domain.numberofcells[1]
Nz = F.domain.numberofcells[2]
dr = F.domain.cellsize[1]
dy = F.domain.cellsize[2]
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x
Fx = F.xvalue
Fy = F.yvalue

# define the vector of cell index
row_index = reshape(G[2:Nr+1,2:Nz+1],Nr*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nr+1,:]
Fw = Fx[1:Nr,:]
Fn = Fy[:,2:Nz+1]
Fs = Fy[:,1:Nz]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# compute the divergence
div_x = (re.*Fe - rw.*Fw)./(dr*rp)
div_y = (Fn - Fs)/dy

# define the RHS Vector
RHSdiv = zeros((Nr+2)*(Nz+2))
RHSdivx = zeros((Nr+2)*(Nz+2))
RHSdivy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nr*Nz)
RHSdivx[row_index] = reshape(div_x,Nr*Nz)
RHSdivy[row_index] = reshape(div_y,Nr*Nz)

(RHSdiv, RHSdivx, RHSdivy)

end



# =============== Divergence 2D Radial Term ============================
function divergenceTermRadial2D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nr = F.domain.numberofcells[1]
Ntheta = F.domain.numberofcells[2]
dr = F.domain.cellsize[1]
dtheta = F.domain.cellsize[2]
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x
Fx = F.xvalue
Fy = F.yvalue

# define the vector of cell index
row_index = reshape(G[2:Nr+1,2:Ntheta+1],Nr*Ntheta) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nr+1,:]
Fw = Fx[1:Nr,:]
Fn = Fy[:,2:Ntheta+1]
Fs = Fy[:,1:Ntheta]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# compute the divergence
div_x = (re.*Fe - rw.*Fw)./(dr*rp)
div_y = (Fn - Fs)./(dtheta*rp)

# define the RHS Vector
RHSdiv = zeros((Nr+2)*(Ntheta+2))
RHSdivx = zeros((Nr+2)*(Ntheta+2))
RHSdivy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nr*Ntheta)
RHSdivx[row_index] = reshape(div_x,Nr*Ntheta)
RHSdivy[row_index] = reshape(div_y,Nr*Ntheta)

(RHSdiv, RHSdivx, RHSdivy)

end

# =============== Divergence 3D Term ============================
function divergenceTerm3D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nx = F.domain.numberofcells[1]
Ny = F.domain.numberofcells[2]
Nz = F.domain.numberofcells[3]
dx = F.domain.cellsize[1]
dy = F.domain.cellsize[2]
dz = F.domain.cellsize[3]
Fx = F.xvalue
Fy = F.yvalue
Fz = F.zvalue

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nx+1,:,:]
Fw = Fx[1:Nx,:,:]
Fn = Fy[:,2:Ny+1,:]
Fs = Fy[:,1:Ny,:]
Ff = Fz[:,:,2:Nz+1]
Fb = Fz[:,:,1:Nz]

# compute the divergence
div_x = (Fe - Fw)/dx
div_y = (Fn - Fs)/dy
div_z = (Ff - Fb)/dz

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)

(RHSdiv, RHSdivx, RHSdivy, RHSdivz)

end



# =============== Divergence 3D Cylindrical Term ============================
function divergenceTermCylindrical3D(F::FaceValue)
# This function calculates the divergence of a field 
# using its face

# extract data from the mesh structure
G = F.domain.numbering
Nx = F.domain.numberofcells[1]
Ny = F.domain.numberofcells[2]
Nz = F.domain.numberofcells[3]
dx = F.domain.cellsize[1]
dy = F.domain.cellsize[2]
dz = F.domain.cellsize[3]
rp = F.domain.cellcenters.x
Fx = F.xvalue
Fy = F.yvalue
Fz = F.zvalue

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the 
# code readability
Fe = Fx[2:Nx+1,:,:]
Fw = Fx[1:Nx,:,:]
Fn = Fy[:,2:Ny+1,:]
Fs = Fy[:,1:Ny,:]
Ff = Fz[:,:,2:Nz+1]
Fb = Fz[:,:,1:Nz]

# compute the divergence
div_x = (Fe - Fw)/dx
div_y = (Fn - Fs)./(dy*rp)
div_z = (Ff - Fb)/dz

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)

(RHSdiv, RHSdivx, RHSdivy, RHSdivz)

end


# ===================== Gradient Term ============================
function gradientTerm(phi)
# calculates the gradient of a variable
# the output is a face variable
d=phi.domain.dimension
if d==1 || d==1.5
  dx = phi.domain.cellsize[1]
  FaceValue(phi.domain,
    (phi.value[2:end]-phi.value[1:end-1])/dx,
    [1.0],
    [1.0])
elseif d==2 || d==2.5
  dx = phi.domain.cellsize[1]
  dy = phi.domain.cellsize[2]
  FaceValue(phi.domain,
    (phi.value[2:end,2:end-1]-phi.value[1:end-1,2:end-1])/dx,
    (phi.value[2:end-1,2:end]-phi.value[2:end-1,1:end-1])/dy,
    [1.0])
elseif d==2.8
  dx = phi.domain.cellsize[1]
  dtheta = phi.domain.cellsize[2]
  Ntheta = phi.domain.numberofcells[2]
  rp = phi.domain.cellcenters.x
  FaceValue(phi.domain,
    (phi.value[2:end,2:end-1]-phi.value[1:end-1,2:end-1])/dx,
    (phi.value[2:end-1,2:end]-phi.value[2:end-1,1:end-1])./(dtheta*rp),
    [1.0])
elseif d==3
  dx = phi.domain.cellsize[1]
  dy = phi.domain.cellsize[2]
  dz = phi.domain.cellsize[3]
  FaceValue(phi.domain,
    (phi.value[2:end,2:end-1,2:end-1]-phi.value[1:end-1,2:end-1,2:end-1])/dx,
    (phi.value[2:end-1,2:end,2:end-1]-phi.value[2:end-1,1:end-1,2:end-1])/dy,
    (phi.value[2:end-1,2:end-1,2:end]-phi.value[2:end-1,2:end-1,1:end-1])/dz)
elseif d==3.2
  dx = phi.domain.cellsize[1]
  dy = phi.domain.cellsize[2]
  dz = phi.domain.cellsize[3]
  Ntheta = phi.domain.numberofcells[2]
  rp = phi.domain.cellcenters.x
  FaceValue(phi.domain,
    (phi.value[2:end,2:end-1,2:end-1]-phi.value[1:end-1,2:end-1,2:end-1])/dx,
    (phi.value[2:end-1,2:end,2:end-1]-phi.value[2:end-1,1:end-1,2:end-1])./(dy*rp),
    (phi.value[2:end-1,2:end-1,2:end]-phi.value[2:end-1,2:end-1,1:end-1])/dz)    
    
end
end