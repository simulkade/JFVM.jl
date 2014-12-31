# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 27 December, 2014
# ===============================

# ====================== 1D CARTESIAN MESH =======================
function createMesh1D(Nx::Int64, Width::Real)
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
		reshape([1:(Nx+2)], Nx+2), 
		[Nx], 
		[dx], 
		CellLocation([1:Nx]*dx.-dx/2,[0.0],[0.0]), 
		FaceLocation([0:Nx]*dx,[0.0],[0.0]),
		[1],
		[1])
end

# ================= 1D CYLINDRICAL MESH ==========================
function createMeshCylindrical1D(Nr::Int64, Radius::Real)
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
		reshape([1:(Nr+2)], Nr+2), 
		[Nr], 
		[dr], 
		CellLocation([1:Nr]*dr.-dr/2,[0.0],[0.0]), 
		FaceLocation([0:Nr]*dr,[0.0],[0.0]),
		[1],
		[1])
end

# ========================= 2D CARTESIAN MESH ==========================
function createMesh2D(Nx::Int64, Ny::Int64, Width::Real, Height::Real)
# builds a uniform 2D mesh:
# Nx is the number of cells in x (horizontal) direction
# Width is the domain length in x direction

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at 
# the boundaries
# cell size is dx
dx = Width/Nx
dy = Height/Ny
G=reshape([1:(Nx+2)*(Ny+2)], Nx+2, Ny+2)
MeshStructure(2, 
	G, 
	[Nx, Ny], 
	[dx, dy], 
	CellLocation([1:Nx]*dx.-dx/2.0,[1:Ny]*dy.-dy/2.0,[0.0]), 
	FaceLocation([0:Nx]*dx,[0:Ny]*dy,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

# ===================== 2D RADIAL MESH ==========================
function createMeshRadial2D(Nr::Int64, Ntheta::Int64, Radius::Real, Angle::Real)
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
end
dr = Radius/Nr
dtheta = Angle/Ntheta
G=reshape([1:(Nr+2)*(Ntheta+2)], Nr+2, Ntheta+2)
MeshStructure(2.8, 
	G, 
	[Nr, Ntheta], 
	[dr, dtheta], 
	CellLocation([1:Nr]*dr.-dr/2.0,[1:Ntheta]*dtheta.-dtheta/2.0,[0.0]), 
	FaceLocation([0:Nr]*dr,[0:Ntheta]*dtheta,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

# ===================== 2D CYLINDRICAL MESH ==========================
function createMeshCylindrical2D(Nr::Int64, Nz::Int64, Radius::Real, Height::Real)
# builds a uniform cylindrical 2D mesh:
# Nr is the number of cells in r (radial) direction
# Radius is the domain length in r direction

# numbering system of cells, like the single index numbering of Matlab
# +2 is added to account for the ghost cells that are added at 
# the boundaries
# cell size is dr, dz
dr = Radius/Nr
dz = Height/Nz
G=reshape([1:(Nr+2)*(Nz+2)], Nr+2, Nz+2)
MeshStructure(2.5, 
	G, 
	[Nr, Nz], 
	[dr, dz], 
	CellLocation([1:Nr]*dr.-dr/2.0,[1:Nz]*dz.-dz/2.0,[0.0]), 
	FaceLocation([0:Nr]*dr,[0:Nz]*dz,[0.0]),
	G[[1,end],[1,end]][:],
	[1])
end

# ========================== 3D CARTESIAN MESH =============================
function createMesh3D(Nx::Int64, Ny::Int64, Nz::Int64, Width::Real, Height::Real, Depth::Real)
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
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2)], Nx+2, Ny+2, Nz+2)
MeshStructure(3, 
	G, 
	[Nx, Ny, Nz], 
	[dx, dy, dz], 
	CellLocation([1:Nx]*dx.-dx/2.0, [1:Ny]*dy.-dy/2.0, [1:Nz]*dz.-dz/2.0), 
	FaceLocation([0:Nx]*dx, [0:Ny]*dy, [0:Nz]*dz),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ny+1, [1, end]][:];
	G[2:Nx+1, [1, end], [1, end]][:]])
end

# ========================== 3D CYLINDRICAL MESH =============================
function createMeshCylindrical3D(Nr::Int64, Ntheta::Int64, Nz::Int64, Radius::Real, Angle::Real, Height::Real)
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
end
dr = Radius/Nr
dtheta = Angle/Ntheta
dz = Height/Nz
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2)], Nr+2, Ntheta+2, Nz+2)
MeshStructure(3.2, 
	G, 
	[Nr, Ntheta, Nz], 
	[dr, dtheta, dz], 
	CellLocation([1:Nr]*dr.-dr/2.0, [1:Ntheta]*dtheta.-dtheta/2.0, [1:Nz]*dz.-dz/2.0), 
	FaceLocation([0:Nr]*dr, [0:Ntheta]*dtheta, [0:Nz]*dz),
	G[[1,end],[1,end],[1,end]][:],
	[G[[1, end], [1, end], 2:Nz+1][:];
	G[[1, end], 2:Ntheta+1, [1, end]][:];
	G[2:Nr+1, [1, end], [1, end]][:]])
end
