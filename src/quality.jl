module quality

using LinearAlgebra

export rankdata, qvalues

function rankdata( values_in::Array{Float64,1} , method::AbstractString="ordinal" )
    if method == "ordinal"
        return( sortperm(values_in) )
    else
        return( values_in )
    end	
end

function fractional_ranks( values_in::Array{Float64,1} , method::AbstractString="ordinal" )
    if method == "ordinal"
        return( (rankdata(values_in,method).-0.5)./length(values_in) )
    else
        return( values_in )
    end	
end

function qvalues( p_values_in::Array{Float64,1} ,
                  pi0::Float64 = 1.0 ,
		  method::AbstractString="ordinal" )
    p_s = p_values_in
    qs_ = []
    m = length(p_s)
    ps = p_s
    frp_  = fractional_ranks( ps,method )
    ifrp_ = [];for (p,f) in zip(ps,frp_) append!(ifrp_, p<=f ? f : p ) end
    #ifrp_ = [ ( (p<=f)*f + p*(p>f) ) for p,f in zip(ps,frp_) ]
    for ip in 1:length(ps)
        p_ = ps[ ip ] ; f_ = frp_[ip]
        q_ = pi0 * p_ / ifrp_[ip]
        append!(qs_, (q_,p_) )
    end
    return(qs_)
end

end
