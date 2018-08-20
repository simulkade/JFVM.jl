# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ====================== 1D CARTESIAN MESH =======================
"""
m = createMesh1D(Nx::Int, Width::Real)
It creates a uniform mesh on a 1D Cartesian domain.

Inputs:

   + Nx: number of cells in the domain (integer)
	 + Width: width or length of the domain (Real)

Outputs:

   + m: a mesh structure

Usage:
```
Nx=10
Lx=1.0
m=createMesh1D(Nx, Lx)
```
"""
function createMesh1D(Nx::Int, Width::Real)
# builds a uniform 1D mesh:
# Nx is the number of cells in x (horizontal) direction
# Width is the domain length in x direction

# mesh dimension

# cell size is dx
dx = Width/Nx

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
MeshStructure(1,
		[Nx],
		CellSize(dx*ones(Nx+2), [0.0], [0.0]),
		CellLocation([1:Nx;]*dx.-dx/2,[0.0],[0.0]),
		FaceLocation([0:Nx;]*dx,[0.0],[0.0]),
		[1],
		[1])
end

"""
m = createMesh1D(facelocationX::Array{Real,1})
It creates a non/uniform mesh on a 1D Cartesian domain.

Inputs:

   + facelocationX: location of the the cell boundaries on
	 the x-axis.

Outputs:

   + m: a mesh structure

Usage:
```
	x= [0.0, 1.0, 1.4, 2.5, 4.1, 6.0, 10.0]
	m=createMesh1D(x)
```
"""
function createMesh1D(facelocationX::Array{T,1}) where T<:Real
# builds a uniform 1D mesh:
# facelocationX is the location of each cell face
# mesh dimension
# cell size is dx
Nx = length(facelocationX)-1
# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
MeshStructure(1,
		[Nx],
		CellSize([facelocationX[2]-facelocationX[1];
		facelocationX[2:end]-facelocationX[1:end-1];
		facelocationX[end]-facelocationX[end-1]], [0.0], [0.0]),
		CellLocation(0.5*(facelocationX[2:end]+facelocationX[1:end-1]),[0.0],[0.0]),
		FaceLocation(facelocationX,[0.0],[0.0]),
		[1],
		[1])
end



# ================= 1D CYLINDRICAL MESH ==========================
"""
m = createMeshCylindrical1D(Nr::Int, Radius::Real)
It creates a uniform mesh on a 1D Radial domain.

Inputs:

   + Nr: number of cells in the domain (integer)
	 + Radius: Radius of the domain (Real)

Outputs:

   + m: a mesh structure

Usage:

Nr=10
R=1.0
m=createMesh1D(Nr, R)
"""
function createMeshCylindrical1D(Nr::Int, Radius::Real)
# builds a uniform 1D cylindrical mesh:
# Nx is the number of cells in r (radial) direction
# Radius is the domain length in r direction
# mesh dimension
# cell size is dr
dr = Radius/Nr
# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
MeshStructure(1.5,
		[Nr],
		CellSize(dr*ones(Nr+2), [0.0], [0.0]),
		CellLocation([1:Nr;]*dr.-dr/2.0,[0.0],[0.0]),
		FaceLocation([0:Nr;]*dr,[0.0],[0.0]),
		[1],
		[1])
end

"""
m = createMeshCylindrical1D(facelocationR::Array{Real,1})
It creates a non/uniform mesh on a 1D radial domain.

Inputs:

   + facelocationX: location of the the cell boundaries on
	 the r-axis.

Outputs:

   + m: a mesh structure

Usage:

	r= [0.1, 1.0, 1.4, 2.5, 4.1, 6.0, 10.0]
	m=createMeshCylindrical1D(r)
"""
function createMeshCylindrical1D(facelocationR::Array{T,1}) where T<:Real
# builds a uniform 1D cylindrical mesh:
# Nx is the number of cells in r (radial) direction
# Radius is the domain length in r direction
# mesh dimension
# cell size is dr
Nr = length(facelocationR)-1
# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
MeshStructure(1.5,
		[Nr],
		CellSize([facelocationR[2]-facelocationR[1];
		facelocationR[2:end]-facelocationR[1:end-1];
		facelocationR[end]-facelocationR[end-1]], [0.0], [0.0]),
		CellLocation(0.5*(facelocationR[2:end]+facelocationR[1:end-1]),[0.0],[0.0]),
		FaceLocation(facelocationR,[0.0],[0.0]),
		[1],
		[1])
end



# ========================= 2D CARTESIAN MESH ==========================
"""
m = createMesh2D(Nx::Int, Ny::Int, Width::Real, Height::Real)
It creates a uniform mesh on a 2D Cartesian domain.

Inputs:

   + Nx: number of cells in the x-direction (integer)
   + Ny: number of cells in the y-direction (integer)
	 + Width: length of the domain in the x-direction (Real)
   + Height: length of the domain in the y-direction (Real)

Outputs:

   + m: a mesh structure

Usage:

Nx=10
Ny=15
Lx=1.0
Ly=2.5
m=createMesh2D(Nx, Ny, Lx, Ly)
"""
function createMesh2D(Nx::Int, Ny::Int, Width::Real, Height::Real)
# builds a uniform 2D mesh:
# Nx is the number of cells in x (horizontal) direction
# Width is the domain length in x direction

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
# cell size is dx
dx = Width/Nx
dy = Height/Ny
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
MeshStructure(2,
	[Nx, Ny],
	CellSize(dx*ones(Nx+2), dy*ones(Ny+2), [0.0]),
	CellLocation([1:Nx;]*dx.-dx/2.0,[1:Ny;]*dy.-dy/2.0,[0.0]),
	FaceLocation([0:Nx;]*dx,[0:Ny;]*dy,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

function createMesh2D(facelocationX::Array{T,1}, facelocationY::Array{T,1}) where T<:Real
Nx = length(facelocationX)-1
Ny = length(facelocationY)-1
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
MeshStructure(2,
	[Nx, Ny],
	CellSize([facelocationX[2]-facelocationX[1]; facelocationX[2:end]-facelocationX[1:end-1]; facelocationX[end]-facelocationX[end-1]],
	[facelocationY[2]-facelocationY[1]; facelocationY[2:end]-facelocationY[1:end-1]; facelocationY[end]-facelocationY[end-1]],
	[0.0]),
	CellLocation(0.5*(facelocationX[2:end]+facelocationX[1:end-1]), 0.5*(facelocationY[2:end]+facelocationY[1:end-1]), [0.0]),
	FaceLocation(facelocationX, facelocationY, [0.0]),
	G[[1,end],[1,end]][:],
	[1])
end



# ============================== 2D RADIAL MESH =======================================
function createMeshRadial2D(Nr::Int, Ntheta::Int, Radius::Real, Angle::Real)
# builds a uniform radial 2D mesh:
# Nr is the number of cells in r (radial) direction
# Radius is the domain length in r direction
# Ntheta is the number of cells in angular direction (not really a direction)
# Angle is the angular domain size
# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
# cell size is dx, dy
if Angle>2*pi
	Angle = 2*pi
	println("The domain size adjusted to match a maximum of 2*pi.")
end
dr = Radius/Nr
dtheta = Angle/Ntheta
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
MeshStructure(2.8,
	[Nr, Ntheta],
	CellSize(dr*ones(Nr+2), dtheta*ones(Ntheta+2), [0.0]),
	CellLocation([1:Nr;]*dr.-dr/2.0,[1:Ntheta;]*dtheta.-dtheta/2.0,[0.0]),
	FaceLocation([0:Nr;]*dr,[0:Ntheta;]*dtheta,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

function createMeshRadial2D(facelocationR::Array{T,1}, facelocationTheta::Array{T,1}) where T<:Real
if facelocationTheta[end]>2.0*pi
	facelocationTheta = facelocationTheta/facelocationTheta[end]*(2.0*pi)
	println("The domain size adjusted to match a maximum of 2*pi.")
end
Nx = length(facelocationR)-1
Ny = length(facelocationTheta)-1
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
MeshStructure(2.8,
	[Nx, Ny],
	CellSize([facelocationR[2]-facelocationR[1]; facelocationR[2:end]-facelocationR[1:end-1]; facelocationR[end]-facelocationR[end-1]],
	[facelocationTheta[2]-facelocationTheta[1]; facelocationTheta[2:end]-facelocationTheta[1:end-1]; facelocationTheta[end]-facelocationTheta[end-1]],
	[0.0]),
	CellLocation(0.5*(facelocationR[2:end]+facelocationR[1:end-1]), 0.5*(facelocationTheta[2:end]+facelocationTheta[1:end-1]), [0.0]),
	FaceLocation(facelocationR, facelocationTheta, [0.0]),
	G[[1,end],[1,end]][:],
	[1])
end



# ===================== 2D CYLINDRICAL MESH ==========================
function createMeshCylindrical2D(Nr::Int, Nz::Int, Radius::Real, Height::Real)
# builds a uniform cylindrical 2D mesh:
# Nr is the number of cells in r (radial) direction
# Radius is the domain length in r direction

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
# cell size is dr, dz
dr = Radius/Nr
dz = Height/Nz
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
MeshStructure(2.5,
	[Nr, Nz],
	CellSize(dr*ones(Nr+2), dz*ones(Nz+2), [0.0]),
	CellLocation([1:Nr;]*dr.-dr/2.0,[1:Nz;]*dz.-dz/2.0,[0.0]),
	FaceLocation([0:Nr;]*dr,[0:Nz;]*dz,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

function createMeshCylindrical2D(facelocationR::Array{T,1}, facelocationY::Array{T,1}) where T<:Real
Nx = length(facelocationR)-1
Ny = length(facelocationY)-1
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
MeshStructure(2.5,
	[Nx, Ny],
	CellSize([facelocationR[2]-facelocationR[1]; facelocationR[2:end]-facelocationR[1:end-1]; facelocationR[end]-facelocationR[end-1]],
	[facelocationY[2]-facelocationY[1]; facelocationY[2:end]-facelocationY[1:end-1]; facelocationY[end]-facelocationY[end-1]],
	[0.0]),
	CellLocation(0.5*(facelocationR[2:end]+facelocationR[1:end-1]), 0.5*(facelocationY[2:end]+facelocationY[1:end-1]), [0.0]),
	FaceLocation(facelocationR, facelocationY, [0.0]),
	G[[1,end],[1,end]][:],
	[1])
end



# ================================= 3D CARTESIAN MESH ========================================
function createMesh3D(Nx::Int, Ny::Int, Nz::Int, Width::Real, Height::Real, Depth::Real)
# builds a uniform 3D mesh:
# Nx is the number of cells in x direction
# Ny is the number of cells in y direction
# Nz is the number of cells in z direction
# Width is the domain length in x direction
# Height is the domain length in y direction
# Depth is the domain length in z direction
# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
dx = Width/Nx
dy = Height/Ny
dz = Depth/Nz
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
MeshStructure(3,
	[Nx, Ny, Nz],
	CellSize(dx*ones(Nx+2), dy*ones(Ny+2), dz*ones(Nz+2)),
	CellLocation([1:Nx;]*dx.-dx/2.0, [1:Ny;]*dy.-dy/2.0, [1:Nz;]*dz.-dz/2.0),
	FaceLocation([0:Nx;]*dx, [0:Ny;]*dy, [0:Nz;]*dz),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ny+1, [1, end]][:];
	G[2:Nx+1, [1, end], [1, end]][:]])
end

function createMesh3D(facelocationX::Array{T,1}, facelocationY::Array{T,1}, facelocationZ::Array{T,1}) where T<:Real
Nx = length(facelocationX)-1
Ny = length(facelocationY)-1
Nz = length(facelocationZ)-1
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
MeshStructure(3,
	[Nx, Ny, Nz],
	CellSize([facelocationX[2]-facelocationX[1]; facelocationX[2:end]-facelocationX[1:end-1]; facelocationX[end]-facelocationX[end-1]],
	  [facelocationY[2]-facelocationY[1]; facelocationY[2:end]-facelocationY[1:end-1]; facelocationY[end]-facelocationY[end-1]],
	  [facelocationZ[2]-facelocationZ[1]; facelocationZ[2:end]-facelocationZ[1:end-1]; facelocationZ[end]-facelocationZ[end-1]]),
	CellLocation(0.5*(facelocationX[2:end]+facelocationX[1:end-1]),
	  0.5*(facelocationY[2:end]+facelocationY[1:end-1]),
	  0.5*(facelocationZ[2:end]+facelocationZ[1:end-1])),
	FaceLocation(facelocationX, facelocationY, facelocationZ),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ny+1, [1, end]][:];
	G[2:Nx+1, [1, end], [1, end]][:]])
end



# ======================================= 3D CYLINDRICAL MESH ==============================================
function createMeshCylindrical3D(Nr::Int, Ntheta::Int, Nz::Int, Radius::Real, Angle::Real, Height::Real)
# builds a uniform 3D mesh:
# Nr is the number of cells in x direction
# Ntheta is the number of cells in y direction
# Nz is the number of cells in z direction
# Radius is the domain length in x direction
# Angle is the domain length in y direction
# Depth is the domain length in z direction

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at
# the boundaries
if Angle>2*pi
	Angle = 2*pi
	println("The domain size adjusted to match a maximum of 2*pi.")
end
dr = Radius/Nr
dtheta = Angle/Ntheta
dz = Height/Nz
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
MeshStructure(3.2,
	[Nr, Ntheta, Nz],
	CellSize(dr*ones(Nr+2), dtheta*ones(Ntheta+2), dz*ones(Nz+2)),
	CellLocation([1:Nr;]*dr.-dr/2.0, [1:Ntheta;]*dtheta.-dtheta/2.0, [1:Nz;]*dz.-dz/2.0),
	FaceLocation([0:Nr;]*dr, [0:Ntheta;]*dtheta, [0:Nz;]*dz),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ntheta+1, [1, end]][:];
	G[2:Nr+1, [1, end], [1, end]][:]])
end

function createMeshCylindrical3D(facelocationR::Array{T,1}, facelocationTheta::Array{T,1}, facelocationZ::Array{T,1}) where T<:Real
if facelocationTheta[end]>2*pi
	facelocationTheta = facelocationTheta/facelocationTheta[end]*2.0*pi
	println("The domain size adjusted to match a maximum of 2*pi.")
end
Nx = length(facelocationR)-1
Ny = length(facelocationTheta)-1
Nz = length(facelocationZ)-1
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
MeshStructure(3.2,
	[Nx, Ny, Nz],
	CellSize([facelocationR[2]-facelocationR[1]; facelocationR[2:end]-facelocationR[1:end-1]; facelocationR[end]-facelocationR[end-1]],
	  [facelocationTheta[2]-facelocationTheta[1]; facelocationTheta[2:end]-facelocationTheta[1:end-1]; facelocationTheta[end]-facelocationTheta[end-1]],
	  [facelocationZ[2]-facelocationZ[1]; facelocationZ[2:end]-facelocationZ[1:end-1]; facelocationZ[end]-facelocationZ[end-1]]),
	CellLocation(0.5*(facelocationR[2:end]+facelocationR[1:end-1]),
	  0.5*(facelocationTheta[2:end]+facelocationTheta[1:end-1]),
	  0.5*(facelocationZ[2:end]+facelocationZ[1:end-1])),
	FaceLocation(facelocationR, facelocationTheta, facelocationZ),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ny+1, [1, end]][:];
	G[2:Nx+1, [1, end], [1, end]][:]])
end
