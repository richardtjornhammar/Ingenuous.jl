include("clustering.jl")
include("models.jl")

println( models.test_matrix )
println( typeof( models.test_matrix ) )
println( models.test_matrix[CartesianIndex(1, 1)] )
println( clustering.connectivity(models.test_matrix,8.0,true) )
println( clustering.calculate_hierarchy_matrix(models.test_distance_matrix) )