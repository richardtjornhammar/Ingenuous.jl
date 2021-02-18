module clustering

using LinearAlgebra

export function distance_matrix_to_absolute_coordinates

function distance_matrix_to_absolute_coordinates ( D::Array{Float64,2},
	 					 bSquared::Bool=false ,
						 n_dimensions::Int64=2 ):
    # Same function can be found in my impetuous-gfa repo here:
    # https://github.com/richardtjornhammar/impetuous/blob/master/src/impetuous/clustering.py
    # around line 121 (same name)
    if not bSquared
        D = D**2.
    end
    N , M = size(D)
    DIM = n_dimensions
    DIJ = D*0.
    for i in 1:M :
        for j in 1:M :
            DIJ[i,j] = 0.5* (D[i,end]+D[j,end]-D[i,j])
    D = DIJ
    println(D)
    """
    U,S,Vt = np.linalg.svd ( D , full_matrices = True )
    S[DIM:] *= 0.
    Z = np.diag(S**0.5)[:,:DIM]
    xr = np.dot( Z.T,Vt )
    return ( xr )
    """
end

end
