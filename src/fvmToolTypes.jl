type CellLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

type CellSize{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

type FaceLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

type MeshStructure
  dimension::Real
  dims::Array{Int,1}
  cellsize::CellSize
  cellcenters::CellLocation
  facecenters::FaceLocation
  corner::Array{Int,1}
  edge::Array{Int,1}
end

type CellValue{T<:Real}
  domain::MeshStructure
  value::Array{T}
end

type FaceValue{T<:Real}
  domain::MeshStructure
  xvalue::Array{T}
  yvalue::Array{T}
  zvalue::Array{T}
end

type BorderValue{T<:Real}
  a::Array{T}
  b::Array{T}
  c::Array{T}
  periodic::Bool
end

type BoundaryCondition
  domain::MeshStructure
  left::BorderValue
  right::BorderValue
  bottom::BorderValue
  top::BorderValue
  back::BorderValue
  front::BorderValue
end