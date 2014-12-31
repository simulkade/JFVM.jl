# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 28 December, 2014
# ===============================

# ==========================================================
# Last changes:
#  * 2014-12-28: - support for cylindrical 3D and radial 2D
#  * 2014-12-30: - cell boundary functions added;
#                - Debugging and testing
# ==========================================================


# ================================ CREATE BOUNDARY CONDITION ======================= 
function createBC(m::MeshStructure)
# creates a boundary condition structure
d=m.dimension
if d==1 || d==1.5
  BoundaryCondition(m,
      BorderValue([1.0], [0.0], [0.0], false),
      BorderValue([1.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false))
elseif d==2 || d==2.5 || d==2.8
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  BoundaryCondition(m,
      BorderValue(ones(1,Ny), zeros(1,Ny), zeros(1,Ny), false),
      BorderValue(ones(1,Ny), zeros(1,Ny), zeros(1,Ny), false),
      BorderValue(ones(Nx,1), zeros(Nx,1), zeros(Nx,1), false),
      BorderValue(ones(Nx,1), zeros(Nx,1), zeros(Nx,1), false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false))
elseif d==3 || d==3.2
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  Nz = m.numberofcells[3]
  BoundaryCondition(m,
      BorderValue(ones(1,Ny,Nz), zeros(1,Ny,Nz), zeros(1,Ny,Nz), false),
      BorderValue(ones(1,Ny,Nz), zeros(1,Ny,Nz), zeros(1,Ny,Nz), false),
      BorderValue(ones(Nx,1,Nz), zeros(Nx,1,Nz), zeros(Nx,1,Nz), false),
      BorderValue(ones(Nx,1,Nz), zeros(Nx,1,Nz), zeros(Nx,1,Nz), false),
      BorderValue(ones(Nx,Ny,1), zeros(Nx,Ny,1), zeros(Nx,Ny,1), false),
      BorderValue(ones(Nx,Ny,1), zeros(Nx,Ny,1), zeros(Nx,Ny,1), false))
end
end
  
function boundaryConditionTerm(BC::BoundaryCondition)
d = BC.domain.dimension
if d==1 || d==1.5
  boundaryCondition1D(BC)
elseif (d == 2) || (d == 2.5)
  boundaryCondition2D(BC)
elseif (d == 2.8)
  boundaryConditionRadial2D(BC)
elseif (d == 3)
  boundaryCondition3D(BC)
elseif (d == 3.2)
  boundaryConditionCylindrical3D(BC)
end

end

# =============================== BOUNDARY CARTESIAN 1D ============================  
function boundaryCondition1D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for 
# a boundary condition structure
G = BC.domain.numbering
Nx = BC.domain.numberofcells[1]
dx = BC.domain.cellsize[1]

# number of boundary nodes:
nb = 4

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros(Nx+2)

q = 0
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
    # Right boundary
    i = Nx+2
    q=q+1
    ii[q] = G[i]  
    jj[q] = G[i]  
    s[q] = BC.right.b[1]/2.0 + BC.right.a[1]/dx
    q=q+1
    ii[q] = G[i]  
    jj[q] = G[i-1] 
    s[q] = BC.right.b[1]/2.0 - BC.right.a[1]/dx
    BCRHS[G[i]] = BC.right.c[1]

    # Left boundary
    i = 1
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i+1] 
    s[q] = -(BC.left.b[1]/2.0 + BC.left.a[1]/dx)
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i] 
    s[q] = -(BC.left.b[1]/2.0 - BC.left.a[1]/dx)
    BCRHS[G[i]] = -BC.left.c[1]
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
    # Right boundary
    i = Nx+2
    q=q+1
    ii[q] = G[i]  
    jj[q] = G[i]  
    s[q] = 1.0
    q=q+1
    ii[q] = G[i]
    jj[q] = G[2]
    s[q] = -1.0 
    BCRHS[G[i]] = 0.0

    # Left boundary
    i = 1
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i]
    s[q] = 1.0
    q=q+1
    ii[q] = G[i]
    jj[q] = G[Nx+1]
    s[q] = -1.0;
    BCRHS[G[i]] = 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], Nx+2, Nx+2)
(BCMatrix, BCRHS)

end



# =============================== BOUNDARY CARTESIAN 2D ============================
function boundaryCondition2D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for 
# a boundary condition structure
G = BC.domain.numbering
Nx, Ny = tuple(BC.domain.numberofcells...)
dx, dy = tuple(BC.domain.cellsize...)

# number of boundary nodes:
nb = 4*(Nx+Ny+2)

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros((Nx+2)*(Ny+2))

for q = 1:4
  ii[q] = BC.domain.corner[q]
  jj[q] = BC.domain.corner[q]
  s[q] = maximum(BC.top.b/2.0+BC.top.a/dy)
  BCRHS[BC.domain.corner[q]] = 0.0
end
q = 4
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic # non-periodic boundary condition
  # top boundary
  j=Ny+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]  
    s[q] = BC.top.b[i-1]/2.0 + BC.top.a[i-1]/dy
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] = BC.top.b[i-1]/2.0 - BC.top.a[i-1]/dy
    BCRHS[G[i,j]] = BC.top.c[i-1]
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] = -(BC.bottom.b[i-1]/2.0 + BC.bottom.a[i-1]/dy)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.bottom.b[i-1]/2.0 - BC.bottom.a[i-1]/dy)
    BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
  end
elseif BC.top.periodic || BC.bottom.periodic  # periodic boundary condition
  # top boundary
  j=Ny+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,2]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ny+1]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end  
end

if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.right.b[j-1]/2.0 + BC.right.a[j-1]/dx
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] = BC.right.b[j-1]/2.0 - BC.right.a[j-1]/dx
    BCRHS[G[i,j]] = BC.right.c[j-1]
  end
  # Left boundary
  i = 1
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] = -(BC.left.b[j-1]/2.0 + BC.left.a[j-1]/dx)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.left.b[j-1]/2.0 - BC.left.a[j-1]/dx)
    BCRHS[G[i,j]] = -(BC.left.c[j-1])
  end  
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[2,j]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
  # Left boundary
  i = 1;
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+1,j]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
(BCMatrix, BCRHS)
end

# =============================== BOUNDARY Radial 2D ============================
function boundaryConditionRadial2D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for 
# a boundary condition structure
G = BC.domain.numbering
Nx, Ntheta = tuple(BC.domain.numberofcells...)
dx, dtheta = tuple(BC.domain.cellsize...)
rp = BC.domain.cellcenters.x
# number of boundary nodes:
nb = 4*(Nx+Ntheta+2)

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros((Nx+2)*(Ntheta+2))

for q = 1:4
  ii[q] = BC.domain.corner[q]
  jj[q] = BC.domain.corner[q]
  s[q] = maximum(BC.top.b/2.0+BC.top.a./(dtheta*rp))
  BCRHS[BC.domain.corner[q]] = 0.0
end
q = 4
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic # non-periodic boundary condition
  # top boundary
  j=Ntheta+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]  
    s[q] = BC.top.b[i-1]/2.0 + BC.top.a[i-1]/(dtheta*rp[i-1])
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] = BC.top.b[i-1]/2.0 - BC.top.a[i-1]/(dtheta*rp[i-1])
    BCRHS[G[i,j]] = BC.top.c[i-1]
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] = -(BC.bottom.b[i-1]/2.0 + BC.bottom.a[i-1]/(dtheta*rp[i-1]))
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.bottom.b[i-1]/2.0 - BC.bottom.a[i-1]/(dtheta*rp[i-1]))
    BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
  end
elseif BC.top.periodic || BC.bottom.periodic  # periodic boundary condition
  # top boundary
  j=Ntheta+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,2]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ntheta+1]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end  
end

if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.right.b[j-1]/2.0 + BC.right.a[j-1]/dx
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] = BC.right.b[j-1]/2.0 - BC.right.a[j-1]/dx
    BCRHS[G[i,j]] = BC.right.c[j-1]
  end
  # Left boundary
  i = 1
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] = -(BC.left.b[j-1]/2.0 + BC.left.a[j-1]/dx)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.left.b[j-1]/2.0 - BC.left.a[j-1]/dx)
    BCRHS[G[i,j]] = -(BC.left.c[j-1])
  end  
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[2,j]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
  # Left boundary
  i = 1;
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+1,j]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ntheta+2), (Nx+2)*(Ntheta+2))
(BCMatrix, BCRHS)
end


# ===================================== BOUNDARY 3D =====================================
function boundaryCondition3D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for 
# a boundary condition structure
G = BC.domain.numbering
Nx, Ny, Nz = tuple(BC.domain.numberofcells...)
dx, dy, dz = tuple(BC.domain.cellsize...)

# number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 4*((Nx+1)*(Ny+1)+(Nx+1)*(Nz+1)+(Ny+1)*(Nz+1))

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64, nb)
jj = zeros(Int64, nb)
s = zeros(Float64, nb)

# define the RHS column vector
BCRHS = zeros(Float64, (Nx+2)*(Ny+2)*(Nz+2))

# assign value to the corner nodes (useless cells)
q = 1:8
ii[q] = BC.domain.corner
jj[q] = BC.domain.corner
s[q] = 1.0 
BCRHS[BC.domain.corner] = 0.0

# assign values to the edges (useless cells)
q = q[end]+[1:length(BC.domain.edge)]
ii[q] = BC.domain.edge
jj[q] = BC.domain.edge
s[q] = 1.0 
BCRHS[BC.domain.edge] = 0.0

# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.top.b/2.0 + BC.top.a/dy
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j-1,k] 
    s[q] = BC.top.b/2.0 - BC.top.a/dy
    BCRHS[G[i,j,k]] = BC.top.c 

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j+1,k]
    s[q] = -(BC.bottom.b/2.0 + BC.bottom.a/dy) # consider the reverse direction of normal
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.bottom.b/2.0 - BC.bottom.a/dy) # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.bottom.c);
elseif BC.top.periodic || BC.bottom.periodic # periodic
    # top boundary
    j=Ny+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,j,k]
    s[q] = 1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,2,k] 
    s[q] = -1
    BCRHS[G[i,j,k]] = 0.0 

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k] 
    s[q] = 1.0
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,Ny+1,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

if !BC.right.periodic && !BC.left.periodic
    # Right boundary
    i=Nx+2
    j=2:Ny+1
    k=2:Nz+1
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.right.b/2.0 + BC.right.a/dx
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i-1,j,k]
    s[q] = BC.right.b/2.0 - BC.right.a/dx
    BCRHS[G[i,j,k]] = BC.right.c

    # Left boundary
    i = 1
    j=2:Ny+1
    k=2:Nz+1
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i+1,j,k]
    s[q] = -(BC.left.b/2.0 + BC.left.a/dx)
    # consider the reverse direction of normal
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,j,k] 
    s[q] = -(BC.left.b/2.0 - BC.left.a/dx)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.left.c)
elseif BC.right.periodic || BC.left.periodic # periodic
    # Right boundary
    i=Nx+2
    j=2:Ny+1
    k=2:Nz+1
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]  
    s[q] = 1.0
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[2,j,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0

    # Left boundary
    i = 1
    j=2:Ny+1
    k=2:Nz+1
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]  
    s[q] = 1.0
    q = q[end]+[1:Ny*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+1,j,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

if !BC.front.periodic && !BC.back.periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k+1]
    s[q] = -(BC.back.b/2.0 + BC.back.a/dz)  # consider the reverse direction of normal
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.back.b/2.0 - BC.back.a/dz)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.back.c)

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.front.b/2.0 + BC.front.a/dz
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k-1]
    s[q] = BC.front.b/2.0 - BC.front.a/dz
    BCRHS[G[i,j,k]] = BC.front.c
elseif BC.front.periodic || BC.back.periodic  # periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1.0
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+1]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1.0
    q = q[end]+[1:Nx*Ny]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,2]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
    (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
    
(BCMatrix, BCRHS)
end

# ===================================== BOUNDARY Cylindrical 3D =====================================
function boundaryConditionCylindrical3D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for 
# a boundary condition structure
G = BC.domain.numbering
Nx, Ntheta, Nz = tuple(BC.domain.numberofcells...)
dx, dtheta, dz = tuple(BC.domain.cellsize...)
rp = BC.domain.cellcenters.x
#rp = repmat(BC.domain.cellcenters.x', 1, Nz)

# number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 4*((Nx+1)*(Ntheta+1)+(Nx+1)*(Nz+1)+(Ntheta+1)*(Nz+1))

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64, nb)
jj = zeros(Int64, nb)
s = zeros(Float64, nb)

# define the RHS column vector
BCRHS = zeros(Float64, (Nx+2)*(Ntheta+2)*(Nz+2))

# assign value to the corner nodes (useless cells)
q = 1:8
ii[q] = BC.domain.corner
jj[q] = BC.domain.corner
s[q] = 1.0 
BCRHS[BC.domain.corner] = 0.0

# assign values to the edges (useless cells)
q = q[end]+[1:length(BC.domain.edge)]
ii[q] = BC.domain.edge
jj[q] = BC.domain.edge
s[q] = 1.0 
BCRHS[BC.domain.edge] = 0.0

# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.top.b/2.0 + BC.top.a./(dtheta*rp)
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j-1,k] 
    s[q] = BC.top.b/2.0 - BC.top.a./(dtheta*rp)
    BCRHS[G[i,j,k]] = BC.top.c 

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j+1,k]
    s[q] = -(BC.bottom.b/2.0 + BC.bottom.a./(dtheta*rp)) # consider the reverse direction of normal
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.bottom.b/2.0 - BC.bottom.a./(dtheta*rp)) # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.bottom.c)
elseif BC.top.periodic || BC.bottom.periodic # periodic
    # top boundary
    j=Ntheta+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,j,k]
    s[q] = 1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,2,k] 
    s[q] = -1
    BCRHS[G[i,j,k]] = 0.0 

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k] 
    s[q] = 1.0
    q = q[end]+[1:Nx*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,Ntheta+1,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

if !BC.right.periodic && !BC.left.periodic
    # Right boundary
    i=Nx+2
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.right.b/2.0 + BC.right.a/dx
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i-1,j,k]
    s[q] = BC.right.b/2.0 - BC.right.a/dx
    BCRHS[G[i,j,k]] = BC.right.c

    # Left boundary
    i = 1
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i+1,j,k]
    s[q] = -(BC.left.b/2.0 + BC.left.a/dx)
    # consider the reverse direction of normal
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]  
    jj[q] = G[i,j,k] 
    s[q] = -(BC.left.b/2.0 - BC.left.a/dx)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.left.c)
elseif BC.right.periodic || BC.left.periodic # periodic
    # Right boundary
    i=Nx+2
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]  
    s[q] = 1.0
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[2,j,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0

    # Left boundary
    i = 1
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]  
    s[q] = 1.0
    q = q[end]+[1:Ntheta*Nz]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+1,j,k]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

if !BC.front.periodic && !BC.back.periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k+1]
    s[q] = -(BC.back.b/2.0 + BC.back.a/dz)  # consider the reverse direction of normal
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.back.b/2.0 - BC.back.a/dz)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.back.c)

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.front.b/2.0 + BC.front.a/dz
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k-1]
    s[q] = BC.front.b/2.0 - BC.front.a/dz
    BCRHS[G[i,j,k]] = BC.front.c
elseif BC.front.periodic || BC.back.periodic  # periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1.0
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+1]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1.0
    q = q[end]+[1:Nx*Ntheta]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,2]
    s[q] = -1.0
    BCRHS[G[i,j,k]] = 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
    (Nx+2)*(Ntheta+2)*(Nz+2), (Nx+2)*(Ntheta+2)*(Nz+2))
    
(BCMatrix, BCRHS)
end


# =============================== Cell Boundary ====================================
# Assign values to ghost cells
function cellBoundary!(phi::CellValue, BC::BoundaryCondition)
d=phi.domain.dimension
if d ==1 || d==1.5
  cellBoundary1D!(phi, BC);
elseif d==2 || d==2.5
  cellBoundary2D!(phi, BC)
elseif d==2.8
  cellBoundaryRadial2D!(phi, BC)
elseif d==3
  cellBoundary3D!(phi, BC)
elseif d==3.2
  cellBoundaryCylindrical3D!(phi, BC)
end
end

# =========================== 1D Cartesian ================================
function cellBoundary1D!(phi::CellValue, BC::BoundaryCondition)
dx = phi.domain.cellsize[1]
# boundary condition (a d\phi/dx + b \phi = c, a column vector of [d a])
# a (phi(i)-phi(i-1))/dx + b (phi(i)+phi(i-1))/2 = c
# phi(i) (a/dx+b/2) + phi(i-1) (-a/dx+b/2) = c
# Right boundary, i=m+2
# phi(i) (a/dx+b/2) = c- phi(i-1) (-a/dx+b/2)
# Left boundary, i=2
#  phi(i-1) (-a/dx+b/2) = c - phi(i) (a/dx+b/2)
# define the new phi
if !BC.left.periodic && !BC.right.periodic
    phiBC = [(BC.left.c[1]-phi.value[2]*(BC.left.a[1]/dx+BC.left.b[1]/2.0))/(-BC.left.a[1]/dx+BC.left.b[1]/2.0);
             phi.value[2:end-1];
             (BC.right.c[1]-phi.value[end-1]*(-BC.right.a[1]/dx+BC.right.b[1]/2.0))/(BC.right.a[1]/dx+BC.right.b[1]/2.0)]  
else
    phiBC = [phi.value[end]; phi.value[2:end-1]; phi.value[1]]
end
phi.value[:] = phiBC[:]
phi
end

# =========================== 2D Cartesian ================================
function cellBoundary2D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
G = phi.domain.numbering
Nx = phi.domain.numberofcells[1]
Ny = phi.domain.numberofcells[2]
dx = phi.domain.cellsize[1]
dy = phi.domain.cellsize[2]

# define the output matrix
phiBC = zeros(Nx+2, Ny+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    for i = 2:Nx+1
      phiBC[i,j]= (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/dy+BC.top.b[i-1]/2.0))/(BC.top.a[i-1]/dy+BC.top.b[i-1]/2.0)
    end
    # Bottom boundary
    j=1;
    for i = 2:Nx+1
      phiBC[i,j]= (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/dy+BC.bottom.b[i-1]/2.0))/(-BC.bottom.a[i-1]/dy+BC.bottom.b[i-1]/2.0)
    end
else
    # top boundary
    j=Ny+2
    for i = 2:Nx+1
      phiBC[i,j]= phi.value[i,2]
    end
    # Bottom boundary
    j=1
    for i = 2:Nx+1
      phiBC[i,j]= phi.value[i,end-1]
    end
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nx+2
    for j = 2:Ny+1
      phiBC[i,j]= (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dx+BC.right.b[j-1]/2.0))/(BC.right.a[j-1]/dx+BC.right.b[j-1]/2.0)
    end

    # Left boundary
    i = 1
    for j = 2:Ny+1
      phiBC[i,j]= (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dx+BC.left.b[j-1]/2.0))/(-BC.left.a[j-1]/dx+BC.left.b[j-1]/2)
    end
else
    # Right boundary
    i = Nx+2
    for j = 2:Ny+1
      phiBC[i,j]= phi.value[2,j]
    end
    # Left boundary
    i = 1
    for j = 2:Ny+1
      phiBC[i,j]= phi.value[end-1,j]
    end
end
phi.value[:] = phiBC[:]
phi
end

# ======================================== Cell Boundary 3D ===================================
function cellBoundary3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
G = phi.domain.numbering
Nx = phi.domain.numberofcells[1]
Ny = phi.domain.numberofcells[2]
Nz = phi.domain.numberofcells[3]
dx = phi.domain.cellsize[1]
dy = phi.domain.cellsize[2]
dz = phi.domain.cellsize[3]

# define the output matrix
phiBC = zeros(Nx+2, Ny+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1,k]).*(-BC.top.a/dy+BC.top.b/2.0))./(BC.top.a/dy+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2,k]).*(BC.bottom.a/dy+BC.bottom.b/2.0))./(-BC.bottom.a/dy+BC.bottom.b/2.0)
else
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1,j,k]).*(-BC.right.a/dx+BC.right.b/2.0))./(BC.right.a/dx+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2,j,k]).*(BC.left.a/dx+BC.left.b/2.0))./(-BC.left.a/dx+BC.left.b/2.0)
else
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]    
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1]).*(-BC.front.a/dz+BC.front.b/2.0))./(BC.front.a/dz+BC.front.b/2.0)

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2]).*(BC.back.a/dz+BC.back.b/2.0))./(-BC.back.a/dz+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Radial 2D ===========================
function cellBoundaryRadial2D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
G = phi.domain.numbering
Nr = phi.domain.numberofcells[1]
Ntheta = phi.domain.numberofcells[2]
dr = phi.domain.cellsize[1]
dtheta = phi.domain.cellsize[2]
rp = phi.domain.cellcenters.x'

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/(dtheta*rp[i-1])+BC.top.b[i-1]/2.0))/(BC.top.a[i-1]/(dtheta*rp[i-1])+BC.top.b[i-1]/2.0)
    end
    # Bottom boundary
    j=1;
    for i = 2:Nr+1
      phiBC[i,j]= (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/(dtheta*rp[i-1])+BC.bottom.b[i-1]/2.0))/(-BC.bottom.a[i-1]/(dtheta*rp[i-1])+BC.bottom.b[i-1]/2.0)
    end
else
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,2]
    end
    # Bottom boundary
    j=1
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,end-1]
    end
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dr+BC.right.b[j-1]/2.0))/(BC.right.a[j-1]/dr+BC.right.b[j-1]/2.0)
    end

    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dr+BC.left.b[j-1]/2.0))/(-BC.left.a[j-1]/dr+BC.left.b[j-1]/2.0)
    end
else
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[2,j]
    end
    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[end-1,j]
    end
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Cylindrical 3D ===========================
function cellBoundaryCylindrical3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
G = phi.domain.numbering
Nr = phi.domain.numberofcells[1]
Ntheta = phi.domain.numberofcells[2]
Nz = phi.domain.numberofcells[3]
dr = phi.domain.cellsize[1]
dtheta = phi.domain.cellsize[2]
dz = phi.domain.cellsize[3]
#rp = zeros(Nx,1,Nz)
#rp[:,1,:] = repmat(phi.domain.cellcenters.x, 1, Nz)
rp = phi.domain.cellcenters.x

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1,k]).*(-BC.top.a./(dtheta*rp)+BC.top.b/2.0))./(BC.top.a./(dtheta*rp)+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2,k]).*(BC.bottom.a./(dtheta*rp)+BC.bottom.b/2.0))./(-BC.bottom.a./(dtheta*rp)+BC.bottom.b/2.0)
else
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1,j,k]).*(-BC.right.a/dr+BC.right.b/2.0))./(BC.right.a/dr+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2,j,k]).*(BC.left.a/dr+BC.left.b/2.0))./(-BC.left.a/dr+BC.left.b/2.0)
else
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]    
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1]).*(-BC.front.a/dz+BC.front.b/2.0))./(BC.front.a/dz+BC.front.b/2.0)

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2]).*(BC.back.a/dz+BC.back.b/2.0))./(-BC.back.a/dz+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end
