# By Sangaré Mariam, 30 Septembre 2022

#Notice, we can use this code for several days computation. The parameter u indicates the number of days (l the day index)

using JuMP,CPLEX
function callD(typ::Int,Tli::Number,Wshare)
	include("milpF.jl")
	include("aff.jl")
        include("cle"*string(typ)*".jl")

	open("Coube","w") do io
      	        println(io,"This file contains the latex code to generate the curves\n")
	end
       production=zeros(T,u)
       for t in 1:T for l in 1:u production[t,l]=L*round(sum(prod1[i,t,l] for  i in 1:n),digits=2) end end
       
       open("Sol_MILPInd.txt_type"*string(typ)*"_n"*string(n), "w") do io
		println(io, "\n***************** Start *****************\n prod=$production\n")
       end
       
       GainA= Array{Float64}(undef, (n,u))
       Fbar=zeros(n,T,u)          	# The periodic amount of energy allocated to each member (Incentive)
       FinI=zeros(n,T,u)           	# The Individual injection of prosumer
       ConsInd=zeros(n,T,u)   		# The individual consumption
       InjInd=zeros(n,T,u)
       Fsum=zeros(T,u)             	# Total injection
       Fins=zeros(T,u)               	# Share of Fum that is locally consumed
       TCons=zeros(T,u)           	# Total consumption 
       Frest=zeros(T,u)              	# Total remaining energy
       Need=zeros(T,u)              	# Total energy need
       qP=zeros(T,u)  			# Total amount of energy injected into BSS            
       Obj=0
       iter=0
       enD=0
       limI=0
       enN=Wshare
       Tcpu=0
       Stock=zeros(T,u)           	# Total BSS state of charge
       if typ==2 limI=n+2 else limI=2 end
       ca=ones(T) 
       cv=ones(T) 
       ca=Pac.*ca
       cv=Pvc.*cv 

        while iter<limI && enN>=Wshare
                for i in 1:n
                        # solve the members problem
                        objl,cpu,gap,Fout,Fin,ConsInd[i,:,:],GainA[i,:],C,E,q=membre(i,iter,GainA[i,:],Fbar[i,:,:],Tli,Need,Fsum,production,Frest,cv,ca)
                        Obj+=objl       
                        gap*=100
                        Tcpu+=cpu    
                        # We update:
                        #=
                                The total injection Fsum[t,u]
                                The total amount of injection that is locally consumed Fins
                                The individual locally withdraw values FinI
                                The individual injection InjInd
                                The total consumption Tcons                                
                                The total need                               


                        =#     
                        view(Fsum,:,:).+=view(Fout,:,:)
                        view(Fins,:,:).+=view(Fin,:,:)
                        view(FinI,i,:,:).=view(Fin,:,:)
                        view(InjInd,i,:,:).=view(Fout,:,:)
                        view(TCons,:,:).+=view(ConsInd,i,:,:)   
                        view(Need,:,:).+=view(C,:,:)
                        for t in 1:T 
                                for l in 1:u 
                                        Stock[t,l]+=round(sum(E[b,t,l]/dist for b in 1:B),digits=2) # the total BSS state of charge
                                        for b in 1:B 
                                                if q[b,t,l]>=0 
                                                        qP[t,l]+=round(q[b,t,l],digits=2)           # the total injection in the storage unit (amount injected in BSS)
                                                end 
                                        end 
                                end 
                        end


                        open("Sol_MILPInd.txt_type"*string(typ)*"_n"*string(n), "a") do io
                                println(io, "iteration $iter , individu $i & $objl & $gap")
                        end
                end
                            
                open("Sol_MILPInd.txt_type"*string(typ)*"_n"*string(n), "a") do io
                        println(io, "\nInteration $iter\n The total objective is : $Obj") 
                end
                som=dist*sum(Fins[t,l] for t in 1:T,l in 1:u)                              # calcule the total amount of energy locally consumed over the horizon
                enN=dist*sum(Fsum[t,l] for t in 1:T,l in 1:u)                              # calcule the total amount of energy injected by prosumers over the horizon

                if iter<2 &&  enN>=Wshare
                        Fbar= KeyF(Fsum,TCons,ConsInd,InjInd,Invest)                        # calculate the incentive to demand response
                else 
                        Frest[:,:].=round.(Fsum[:,:].-Fins[:,:],digits=2)                   # Calculate the remaining energy
                        Fbar=Maj(FinI,Fbar,Frest,InjInd)
                end 

                enD=dist*sum(Fbar[i,t,l] for i in 1:n,t in 1:T, l in 1:u)
                perte=enN-som                                                               # calculate the total amount of energy injected into the main grid
                open("Sol_MILPInd.txt_type"*string(typ)*"_n"*string(n), "a") do io
                        println(io, "\nThe available energy in kWh is $enN \nEnergy loss is $perte in kWh\nInj=",round.(Fsum,digits=2),"\nTCons=",round.(TCons,digits=2),"\nFrest=$Frest\nNeed=",round.(Need,digits=2),"\nStock=",round.(Stock,digits=2),"\n Fins=",round.(Fins,digits=2))
                end

                if iter==limI-1 
                        Frest[:,:].=Fsum[:,:]-Fins[:,:]						# Total periodic injection minus the the total lical comsumption
                        affiche(Fsum,Stock,Frest,Need,production,TCons,qP,Fins,typ) # afficher les courbes
                end
                Obj=0
                TCons[:,:].=0
                Fsum[:,:].=0
                Stock[:,:].=0
                Fins[:,:].=0
                Need[:,:].=0
                qP[:,:].=0
                iter+=1
        end                                                                			 #End While

        open("Sol_MILPInd.txt_type"*string(typ)*"_n"*string(n), "a") do io
                println(io, "******************************** Distributed is finished, the total CPU is $Tcpu ******************************************************\n")
        end

end

function callF()
                                                # typ determine the type of reparatition keys if type==1 SWEEN, 
        Tli=50.0
        Nbcle=5                                 # Number of key calculation methods
        Wshare=0.0                           	# Is the threshold of injection beyond which we share the energy. If there is not enough energy there is no need to share.
        for n in [7]
                @info "************************************************ Main : $n member instance *****************************************"
                include("Data/n"*string(n)*"_t24.txt")
                #include("Data/a.jl")
                for typ in 1:Nbcle
                        callD(typ,Tli,Wshare)
                end
        end
end


 callF()
  




