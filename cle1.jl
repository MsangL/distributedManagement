
# Calcul clé de répartition selon le modèle existant
function KeyF(Fsum::Matrix{Float64},Tcons::Matrix{Float64},ConsInd::Array{Float64},InjInd::Array{Float64},Invest::Matrix{Float64})
      keyI=zeros(n,T,u)
     indexI=zeros(n)
     uu=0
     for l in 1:u
            for t in 1:T
                  if Fsum[t,l]>0                                                   # If the total injection is positive                                                  
                        if Fsum[t,l]>Tcons[t,l]                                    # If the total injection is greater than the total consumption
                              for i in 1:n
                                    keyI[i,t,l]=ConsInd[i,t,l]                     # Each member receive keyI[i,t,l] equal to her need
                              end
                        else                                                       # if total production is less then total consumtion
                              Cmin=MinimumCons(ConsInd,t,l)                        # find minimum consumption 
                              if n*Cmin<=Fsum[t,l]                                 # If each member can receive the minimum amount Cmin
                                    
                                    for i in 1:n
                                          indexI[i]=i                              # indexI is an auxillary vector which indicates if a member's need is satisfied. A member i is satisfied if index_i=0
                                    end
                                      
                                    uu=n                                           # the initial unsatisfied members
                                    residu=Fsum[t,l]
                                    while uu*Cmin<=residu                          # while we can allocate to each unsatisfied member, 
                                          numb=0
                                          for i in 1:n
                                                if indexI[i]!=0
                                                      keyI[i,t,l]+=Cmin            # update each member's amount with Cmin           
                                                end
                                          end
                                          for i in 1:n 
                                                if ConsInd[i,t,l]==Cmin             # The members whose need was Cmin are satisfied, we determine the number of satisfied members, we update their consumption with a large enough value to ensure their Cmin is always < to their consumption, we set their indexI equals to 0
                                                     numb+=1
                                                     indexI[i]=0
                                                     ConsInd[i,t,l]=+10000
                                                end
                                          end
                                          residu-=uu*Cmin                           # we calculate the ramining energy
                                          uu-=numb                                  # we update the number of unsatisfied members
                                          Cmin=MinimumCons(ConsInd,t,l)             #  we determine Cmin for the next iteration
                                    end                                                                      # end Do while                                           
                                    for i in 1:n
                                          if indexI[i]!=0
                                                keyI[i,t,l]+=(residu/uu)            # After leaving the while, we uniformly share the remaining energy between the remaing members 
                                          end
                                          #keyI[i,t,l]
                                    end
                              else                                                                                               #       n*Cmin>Fsum                     
                                    for i in 1:n
                                          keyI[i,t,l]=Fsum[t,l]/n                   # if we cannot allocate we initial Cmin to each member, we share uniformly
                                    end
                              end
                        end
                  else                                                                                  # Fsum=0
                        for i in 1:n
                              keyI[i,t,l]=0                                         # if there is no energy to share
                        end
                  end
            end                                                                          # end for t
     end                                                                                 # end for l
      return keyI
end


function MinimumCons(ConsInd::Array{Float64},t::Int,l::Int)
      Cmin=+10000
      for i in 1:n
            if ConsInd[i,t,l]<=Cmin
                    Cmin=ConsInd[i,t,l]
            end
      end
      return Cmin
end


