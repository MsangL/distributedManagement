infi=10000

#refer to cle1.jl for more details
function KeyF(Fsum::Matrix{Float64},Tcons::Matrix{Float64},ConsInd::Array{Float64},InjInd::Array{Float64},Invest::Matrix{Float64})
      cand= ones(n,T,u)                  # this will indicate if a member can receive energy
      Tcand= n.*ones(T,u)          # this is the number of member  in need 
      keyI= zeros(n,T,u)
      indexI= Array{Float64}(undef, (n))

      for i in 1:n
            for t in 1:T
                  for l in 1:u
                        if InjInd[i,t,l]>0
                              cand[i,t,l]=0
                              Tcand[t,l]-=1
                              ConsInd[i,t,l]=infi
                        end
                  end
            end
      end

     uu=0
     for l in 1:u
            for t in 1:T
                  
                if Fsum[t,l]>0                                                                                                                           # si les injections sont possibles
                         if Fsum[t,l]>sum(ConsInd[i,t,l] for i in 1:n if cand[i,t,l]>0)                                    # if total production is more then total consumption
                              for i in 1:n
                                    keyI[i,t,l]=ConsInd[i,t,l]*cand[i,t,l]
                              end
                         else                                                                        # if total production is less then total consumtion
                               Cmin=MinimumCons(ConsInd,t,l)           # find minimum consumption
                               if Tcand[t,l]*Cmin<=Fsum[t,l]
                                      for i in 1:n
                                                 if cand[i,t,l]!=0
                                                       indexI[i]=i
                                                 else
                                                       indexI[i]=0
                                                 end
                                      end
                                      uu=Tcand[t,l]
                                      residu=Fsum[t,l]
                                      while uu*Cmin<=residu 
                                          numb=0
                                          for i in 1:n
                                                if indexI[i]!=0 
                                                      keyI[i,t,l]+=Cmin
                                                end
                                          end
                                          for i in 1:n 
                                                if ConsInd[i,t,l]==Cmin 
                                                     numb+=1
                                                     indexI[i]=0
                                                     ConsInd[i,t,l]=+infi
                                                end
                                          end
                                          residu-=uu*Cmin
                                          uu-=numb
                                          Cmin=MinimumCons(ConsInd,t,l)           
                                     end                                                                      # end Do while                                           
                                     for i in 1:n
                                            if indexI[i]!=0 
                                                  keyI[i,t,l]+=(residu/uu)
                                            end
                                     end
                               else                                                                                               #       n*Cmin>Fsum                     
                                     for i in 1:n
                                           keyI[i,t,l]=(Fsum[t,l]*cand[i,t,l])/Tcand[t,l]
                                     end
                               end
                               
                         end
                   else                                                                                  # Fsum=0
                          for i in 1:n
                                keyI[i,t,l]=0
                          end
                 end                                                                      # end Fsum>0
                   
            end                                                                          # end for t
     end                                                                                 # end for l
       return keyI
end


function MinimumCons(ConsInd::Array{Float64},t::Int,l::Int)
      Cmin=+infi
      for i in 1:n
            if ConsInd[i,t,l]<=Cmin
                    Cmin=ConsInd[i,t,l]
            end
      end
      return Cmin
end

