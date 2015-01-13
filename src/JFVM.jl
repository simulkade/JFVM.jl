module JFVM

using PyPlot, PyCall
pygui_start(:wx)
@pyimport mayavi.mlab as mayavis

export MeshStructure, BoundaryCondition, CellValue, FaceValue,
       arithmeticMean, geometricMean, harmonicMean, upwindMean,
       tvdMean, createBC, boundaryConditionTerm, cellBoundary,
       divergenceTerm, gradientTerm, convectionUpwindTerm, 
       convectionTerm, convectionTvdTerm, diffusionTerm, createCellVariable,
       createFaceVariable, copyCell, fluxLimiter, createMesh1D,
       createMesh2D, createMesh3D, createMeshRadial2D, createMeshCylindrical2D,
       createMeshCylindrical3D, createMeshCylindrical1D, solveLinearPDE, visualizeCells,
       linearSourceTerm, constantSourceTerm, transientTerm, solveMUMPSLinearPDE

include("fvmToolTypes.jl")
include("meshstructure.jl")
include("boundarycondition.jl")
include("domainVariables.jl")
include("diffusionterms.jl")
include("transientTerms.jl")
include("domainOperators.jl")
include("convectionTerms.jl")
include("averagingTerms.jl")
include("calculusTerms.jl")
include("sourceTerms.jl")
include("solveVisualizePDE.jl")
include("juliaFVMtools.jl")

end # module
