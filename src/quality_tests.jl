include("quality.jl")
include("models.jl")

using Random
rvec = rand(10)
println( quality.rankdata( rvec ) )
println( quality.fractional_ranks( rvec ) )
println( quality.qvalues(models.pvalue_list) )