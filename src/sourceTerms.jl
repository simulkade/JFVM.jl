# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2015-01-10 changed numberofcells to dims
# ================================================================

# ======================= Linear source term ========================
function linearSourceTerm(betta0::CellValue)
m = betta0.domain
d = m.dimension
if (d ==1) || (d==1.5)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  b = betta0.value[2:end-1]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  b = betta0.value[2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif (d == 3) || (d==3.2)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  b = betta0.value[2:end-1,2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


function linearSourceTerm(m::MeshStructure, betta0::Real)
linearSourceTerm(CellValue(m, betta0))
end


function linearSourceTerm{T<:Real}(m::MeshStructure, betta0::Array{T})
d = m.dimension
if (d ==1) || (d==1.5)
  Nx = m.dims
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif (d == 3) || (d==3.2)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


# ================================== constant source term ================================
function constantSourceTerm(phi0::CellValue)
m = phi0.domain
d = m.dimension
if (d ==1) || (d==1.5)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] = reshape(phi0.value[2:end-1],Nx)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1],Nx*Ny)
elseif (d == 3) || (d==3.2)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1,2:end-1],Nx*Ny*Nz)
end
RHS
end


function constantSourceTerm(m::MeshStructure, phi0::Real)
constantSourceTerm(CellValue(m, phi0))
end


function constantSourceTerm{T<:Real}(m::MeshStructure, phi0::Array{T})
d = m.dimension
if (d ==1) || (d==1.5)
  Nx = m.dims
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] = reshape(phi0,Nx)
elseif (d == 2) || (d == 2.5) || (d==2.8)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] = reshape(phi0,Nx*Ny)
elseif (d == 3) || (d==3.2)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] = reshape(phi0,Nx*Ny*Nz)
end
RHS
end
