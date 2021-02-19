#import Distances

include("projection.jl")
include("models.jl")

B = vcat(models.test_model,models.test_model.+9);B[:,3]*=0.;B
#D = Distances.pairwise(Distances.Euclidean(),transpose(B),dims=2)
#sum( Distances.pairwise(Distances.Euclidean(),transpose( Ingenuous.projection.distance_matrix_to_absolute_coordinates(D,false) ))-D )

A = projection.distance_matrix_to_absolute_coordinates( models.test_distance_matrix,false )

println( B )
println( )
println( A )
