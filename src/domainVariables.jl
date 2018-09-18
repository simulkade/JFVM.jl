# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 30 December, 2014
# ===============================

# =====================================================================
# 2014-12-30: - cell variables can be created including bounddary cells
# =====================================================================

# ============================================================

"""
creates a cell variable over the domain m and assigns value phi0 to it
function createCellVariable(m::MeshStructure, phi0::Real)
"""
function createCellVariable(m::MeshStructure, phi0::Real)
# creates a cell variable and assigns value phi0 to it
CellValue(m, phi0*ones(tuple(m.dims.+2...)))
end

function createCellVariable(phi::CellValue)
  # similar to copyCell
  CellValue(phi.domain, Base.copy(phi.value))
end

# ============================================================
function createCellVariable(m::MeshStructure, phi0::Array{T}) where T<:Real
# creates a cell variable and assigns value phi0 to it
if prod(m.dims.+2)==length(phi0)
  CellValue(m, phi0)
elseif prod(m.dims)==length(phi0)
  d=m.dimension
  phival = zeros(tuple(m.dims.+2...))
  BC = createBC(m) # Neumann boundaries
  if d==1 || d==1.5
    phival[2:end-1] = phi0
  elseif d==2 d==2 || d==2.5 || d==2.8
    phival[2:end-1, 2:end-1] = phi0
  elseif d==3 || d==3.2
    phival[2:end-1,2:end-1,2:end-1] = phi0
  end
  phi = CellValue(m, phival)
  cellBoundary!(phi, BC)
else
  error("JFVM: Matrix must be the same size as the domain.")
end
end

# ==============================================================
function createCellVariable(m::MeshStructure, phi0::Array{T}, BC::BoundaryCondition) where T<:Real
# creates a cell variable and assigns value phi0 to it
if prod(m.dims.+2)==length(phi0)
  error("JFVM: Matrix must be the same size as the domain.")
elseif prod(m.dims)==length(phi0)
  d=m.dimension
  phival = zeros(tuple(m.dims.+2...))
  if d==1 || d==1.5
    phival[2:end-1] = phi0
  elseif d==2 || d==2.5 || d==2.8
    phival[2:end-1, 2:end-1] = phi0
  elseif d==3 || d==3.2
    phival[2:end-1,2:end-1,2:end-1] = phi0
  end
  phi = CellValue(m, phival)
  cellBoundary!(phi, BC)
else
  error("JFVM: Matrix must be the same size as the domain.")
end
end



# ============================================================
function createCellVariable(m::MeshStructure, phi0::Real, BC::BoundaryCondition)
phi = CellValue(m, phi0*ones(tuple(m.dims.+2...)))
cellBoundary!(phi, BC)
end



# ============================================================
function createFaceVariable(m::MeshStructure, phi0::Array{T,1}) where T<:Real
# creates a face variable based on the mesh structure
d=m.dimension
  if d==1 || d==1.5
    FaceValue(m,
	      ones(m.dims[1]+1)*phi0[1],
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2])*phi0[1],
	      ones(m.dims[1], m.dims[2]+1)*phi0[2],
	      [1.0])
  elseif d==3 || d==3.2
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2], m.dims[3])*phi0[1],
	      ones(m.dims[1], m.dims[2]+1, m.dims[3])*phi0[2],
	      ones(m.dims[1], m.dims[2], m.dims[3]+1)*phi0[3])
  end
end

function createFaceVariable(m::MeshStructure, phi0::Real)
# creates a face variable based on the mesh structure
d=m.dimension
  if d==1 || d==1.5
    FaceValue(m,
	      ones(m.dims[1]+1)*phi0,
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2])*phi0,
	      ones(m.dims[1], m.dims[2]+1)*phi0,
	      [1.0])
  elseif d==3 || d==3.2
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2], m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2]+1, m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2], m.dims[3]+1)*phi0)
  end
end


# ============================================================
function createCellVector(m::MeshStructure, phi0::Array{T,1}) where T<:Real
# creates a cell vector based on the mesh structure
d=m.dimension
  if d==1 || d==1.5
    CellVector(m,
	      ones(m.dims[1])*phi0[1],
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    CellVector(m,
	      ones(m.dims[1], m.dims[2])*phi0[1],
	      ones(m.dims[1], m.dims[2])*phi0[2],
	      [1.0])
  elseif d==3 || d==3.2
    CellVector(m,
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[1],
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[2],
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[3])
  end
end

function createCellVector(m::MeshStructure, phi0::Real)
# creates a cell vector based on the mesh structure
d=m.dimension
  if d==1 || d==1.5
    CellVector(m,
	      ones(m.dims[1])*phi0,
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    CellVector(m,
	      ones(m.dims[1], m.dims[2])*phi0,
	      ones(m.dims[1], m.dims[2])*phi0,
	      [1.0])
  elseif d==3 || d==3.2
    CellVector(m,
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0)
  end
end

# ================== copy function for cell variables ===================
function copyCell(phi::CellValue)
  CellValue(phi.domain, Base.copy(phi.value))
end
