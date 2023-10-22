
# Calcul clé de répartition au prorata de la consommation
function KeyF(Fsum::Matrix{Float64},Tcons::Matrix{Float64},ConsInd::Array{Float64},InjInd::Array{Float64},Invest::Matrix{Float64})
    cand= ones(n,T,u)                  # this will indicate if a member can receive energy
     keyI= zeros(n,T,u)
     
     cand[:,:,:].=1
     
     for i in 1:n
         for t in 1:T
             for l in 1:u
                 if Indnj[i,t,l]>0
                       cand[i,t,l]=0
                 end
             end
         end
     end

	for t in 1:T
		for i in 1:n
			for l in 1:u
				keyI[i,t,l]=(ConsInd[i,t,l]*Fsum[t,l]*cand[i,t,l])/sum(ConsInd[j,t,l] for j in 1:n if cand[j,t,l]>0)
			end
		end
	end
     
       return keyI
end

