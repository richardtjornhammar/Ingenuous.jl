include("clustering.jl")

println( clustering.test_matrix )
println( typeof( clustering.test_matrix ) )
println( clustering.test_matrix[CartesianIndex(1, 1)] )
println( clustering.connectivity(clustering.test_matrix,8.0,true) )