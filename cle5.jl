# Calcul clé de répartition au prorata de l'investissement
function KeyF(Fsum::Matrix{Float64},Tcons::Matrix{Float64},ConsInd::Array{Float64},InjInd::Array{Float64},Invest::Matrix{Float64})
    cand= ones(n,T,u)                  # this will indicate if a member can receive energy
    keyI= zeros(n,T,u)
    
     
     for i in 1:n
         for t in 1:T
             for l in 1:u
                 if InjInd[i,t,l]>0
                       cand[i,t,l]=0
                 end
             end
         end
     end
     
	for t in 1:T
		for i in 1:n
			for l in 1:u
				if Invest[i]>0
					   keyI[i,t,l]=(Invest[i]*Fsum[t,l]*cand[i,t,l])/sum(Invest[j] for j in 1:n if cand[j,t,l]>0)
				else keyI[i,t,l]=0
				end
			end
		end
	end
     
       return keyI
end


