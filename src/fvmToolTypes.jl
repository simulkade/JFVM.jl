type CellLocation
  x::Array{Float64,1}
  y::Array{Float64,1}
  z::Array{Float64,1}
end

type CellSize
  x::Array{Float64,1}
  y::Array{Float64,1}
  z::Array{Float64,1}
end

type FaceLocation
  x::Array{Float64,1}
  y::Array{Float64,1}
  z::Array{Float64,1}
end

type MeshStructure
  dimension::Real
  dims::Array{Int64,1}
  cellsize::CellSize
  cellcenters::CellLocation
  facecenters::FaceLocation
  corner::Array{Int64,1}
  edge::Array{Int64,1}
end

type CellValue
  domain::MeshStructure
  value::Array{Float64}
end

type FaceValue
  domain::MeshStructure
  xvalue::Array{Float64}
  yvalue::Array{Float64}
  zvalue::Array{Float64}
end

type BorderValue
  a::Array{Float64}
  b::Array{Float64}
  c::Array{Float64}
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