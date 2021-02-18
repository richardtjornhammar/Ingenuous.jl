module clustering

using LinearAlgebra

export connectivity,test_matrix

test_matrix = reshape( [ [0.00, 0.10, 0.10, 9.00, 9.00, 9.00] ;
              	[0.10, 0.00, 0.15, 9.00, 9.00, 9.00] ;
        	[0.10, 0.15, 0.00, 9.00, 9.00, 9.00] ;
        	[9.00, 9.00, 9.00, 0.00, 0.10, 0.10] ;
        	[9.10, 9.00, 9.00, 0.10, 0.00, 0.15] ;
        	[9.10, 9.00, 9.00, 0.10, 0.15, 0.00] ] , (6,6) )

function connectivity( B::Array{Float64,2} , val::Float64 , bVerbose::Bool=false )
    """
This is a cutoff based clustering algorithm. The intended use is to supply a distance matrix and a cutoff value (then becomes symmetric positive definite).  For a small distance cutoff, you should see all the parts of the system and for a large distance cutoff, you should see the entire system. It has been employed for statistical analysis work as well as the original application where it was employed to segment molecular systems.
        #    
	# JULIA ADAPTATION OF MY C++ CODE THAT CAN BE FOUND IN
	# https://github.com/richardtjornhammar/RichTools/blob/master/src/cluster.cc
	# AROUND LINE 2277
	# FOR A DESCRIPTION READ PAGE 30 (16 INTERNAL NUMBERING) of:
	# https://kth.diva-portal.org/smash/get/diva2:748464/FULLTEXT01.pdf
	#
    """
    nr_sq,mr_sq = size(B)
    if nr_sq != mr_sq 
        println( "ERROR" )
        return(-1)
    end
    N = mr_sq
    res,nvisi,s,NN,ndx,C = [],[],[],[],[],0
    res = append!(res,[0])
    for i in 1:N
        append!(nvisi,[i])
        append!(res,[0]); append!(res,[0])
        append!(ndx,[i])
    end
    if bVerbose
        println(res," ",ndx," ",nvisi)
        println(nvisi[end])
    end
    while length(ndx)>0
        i = pop!(ndx)
        NN = []
        if nvisi[i]>0
            C-=1
            for j in 1:N
                if B[j,i]<=val
                    append!(NN,[j])
                end
            end
            while length(NN)>0
                k = pop!(NN)
                nvisi[k] = C
                for j in 1:N
                    if B[j,k]<=val
                        for q in 1:N
                            if nvisi[q]==j
                                append!(NN,[q])
                            end
                        end
                    end
                end
            end
        end
    end
    if bVerbose
        println("INFO ", C,"\ncluster data:")
        Nc = []
        for i in 1:-1*C
            append!(Nc,[0])
        end
    end
    for q in 1:N
        res[q*2]=q
        res[q*2-1]=nvisi[q]-C+1
        if bVerbose
	    Nc[res[q*2-1]]+=1
            println(" ",res[q*2-1],"\t",res[2*q])
        end
    end
    if bVerbose
        for i in 1:-1*C
            println("Cluster ",i," has ", Nc[i]," elements")
        end
    end
    pop!(res)
    ret = transpose(reshape(res,(2,N)))
    if bVerbose
        println( ret, size(ret) )
    end
    return ( ret )
end

end
