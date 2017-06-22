module JFVM

using PyPlot
# using PyCall
# I prefer not to use the following command for the issues that it has on windows machines
# pygui_start(:wx)
# mayavis=0
# try
#   @pyimport mayavi.mlab as m
#   mayavis=m
# catch
#   warn("Mayavi is not installed or could not be imported.")
# end

import Base: +, -, *, /
export MeshStructure, BoundaryCondition, CellValue, FaceValue, CellVector,
       arithmeticMean, geometricMean, harmonicMean, upwindMean, linearMean,
       tvdMean, createBC, boundaryConditionTerm, cellBoundary, solvePDE,
       divergenceTerm, gradientTerm, convectionUpwindTerm, createCellVector,
       convectionTerm, convectionTvdTerm, diffusionTerm, createCellVariable,
       createFaceVariable, copyCell, fluxLimiter, createMesh1D,
       createMesh2D, createMesh3D, createMeshRadial2D, createMeshCylindrical2D,
       createMeshCylindrical3D, createMeshCylindrical1D, solveLinearPDE,
       visualizeCells, linearSourceTerm, constantSourceTerm, transientTerm,
       solveMUMPSLinearPDE, faceEval, cellEval, permfieldlogrndg, permfieldlogrnde,
       plot, imshow, xlabel, ylabel, figure, legend, pcolor, contour, colorbar,
       visualizeCellVectors, JFVM_test, solveExplicitPDE, reshapeCell,
       cellVolume, reshapeInternalCell, internalCells, domainInt, convectionTvdRHS

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
include("JFVMtools.jl")
include("jfvm_test.jl")

end # module
