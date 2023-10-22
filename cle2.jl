
function KeyF(Fsum::Matrix{Float64},Tcons::Matrix{Float64},ConsInd::Array{Float64},InjInd::Array{Float64},Invest::Matrix{Float64})
    Fbar= zeros(n,T,u)
    cand= ones(n,T,u)                                     # this will indicate if a member can receive energy
    Tcand=n.*ones(T,u)                                    # this is the number of member  in need 

     
     for i in 1:n
         for t in 1:T
             for l in 1:u
                 if InjInd[i,t,l]>0
                    cand[i,t,l]=0
                    Tcand[t,l]-=1
                 end
             end
         end
     end
     
    for i in 1:n
        for t in 1:T
            for l in 1:u
                if Tcand[t,l]>0
                   	Fbar[i,t,l]=(Fsum[t,l]*cand[i,t,l])/Tcand[t,l]
                end
            end
        end
    end  
      return Fbar
end

               

function Maj(Fin::Array{Float64},Fbar::Array{Float64},reste::Matrix{Float64},inj::Array{Float64})
    cand= ones(n,T,u)                   # this will indicate if a member can receive energy
    Tcand=n.*ones(T,u)                  # this is the number of member  in need 
    Fbar1=zeros(n,T,u)

    reste=zeros(T,u)
    Fbar1[:,:,:].=Fbar[:,:,:]

    for i in 1:n
        for t in 1:T
            for l in 1:u
                if Fbar[i,t,l]-Fin[i,t,l]>0 || inj[i,t,l]>0
                    cand[i,t,l]=0
                    Tcand[t,l]-=1
                    Fbar1[i,t,l]=Fin[i,t,l]
                    reste[t,l]=reste[t,l]+(Fbar[i,t,l]-Fin[i,t,l])
                end
            end
        end

	end         
                     
    for i in 1:n
        for t in 1:T
            for l in 1:u
                if cand[i,t,l]!=0
                    Fbar1[i,t,l]+=(reste[t,l]*cand[i,t,l])/Tcand[t,l]
                end
            end
        end
    end
    return Fbar1
end
