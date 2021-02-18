module projection

using LinearAlgebra

export distance_matrix_to_absolute_coordinates

function distance_matrix_to_absolute_coordinates( D::Array{Float64,2} , bSquared::Bool=false , n_dimensions::Int64=2 )
    """
    Same function can be found in my impetuous-gfa repo here:
    https://github.com/richardtjornhammar/impetuous/blob/master/src/impetuous/clustering.py
    around line 121 (same name)
    """
    if !bSquared
        D = D.^2.
    end
    N,M = size(D)
    DIM = n_dimensions
    DIJ = D.*0.
    for i in 1:M
        for j in 1:M
            DIJ[i,j] = 0.5* (D[i,end]+D[j,end]-D[i,j])
	end
    end
    D = DIJ
    U,S,Vt = svd( D )
    S[DIM+1:end] *= 0.
    Z  = diagm(S.^0.5)[:,1:DIM]
    xr = Vt*Z
    return ( xr )
end

end
