# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 30 December, 2014
# ===============================


# ======================= Linear source term ========================
function linearSourceTerm(betta0::CellValue)
m = betta0.domain
d = m.dimension
G = m.numbering
if (d ==1) || (d==1.5)
  Nx = m.numberofcells[1]
  b = betta0.value[2:end-1]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  b = betta0.value[2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif (d == 3) || (d==3.2)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  Nz = m.numberofcells[3]
  b = betta0.value[2:end-1,2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


function linearSourceTerm(m::MeshStructure, betta0::Real)
linearSourceTerm(CellValue(m, betta0))
end


function linearSourceTerm(m::MeshStructure, betta0::Array{Real})
d = m.dimension
G = m.numbering
if (d ==1) || (d==1.5)
  Nx = m.numberofcells
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif (d == 3) || (d==3.2)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  Nz = m.numberofcells[3]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


# ================================== constant source term ================================
function constantSourceTerm(phi0::CellValue)
m = phi0.domain
d = m.dimension
G = m.numbering
if (d ==1) || (d==1.5)
  Nx = m.numberofcells[1]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] = reshape(phi0.value[2:end-1],Nx)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1],Nx*Ny)
elseif (d == 3) || (d==3.2)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  Nz = m.numberofcells[3]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1,2:end-1],Nx*Ny*Nz)
end
RHS
end


function constantSourceTerm(m::MeshStructure, phi0::Real)
constantSourceTerm(CellValue(m, phi0))
end


function constantSourceTerm(m::MeshStructure, phi0::Array{Real})
d = m.dimension
G = m.numbering
if (d ==1) || (d==1.5)
  Nx = m.numberofcells
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] = reshape(phi0,Nx)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] = reshape(phi0,Nx*Ny)
elseif (d == 3) || (d==3.2)
  Nx = m.numberofcells[1]
  Ny = m.numberofcells[2]
  Nz = m.numberofcells[3]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] = reshape(phi0,Nx*Ny*Nz)
end
RHS
end
