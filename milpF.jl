### By SANGARE Mariam ###

function membre(i::Int,iter::Int,GainA,Fbar::Matrix{Float64},Tli::Number,Need::Matrix{Float64},Fsum::Matrix{Float64},prod::Matrix{Float64},Frest::Matrix{Float64},Pvc::Vector{Float64},Pac::Vector{Float64})
      mod=Model(CPLEX.Optimizer)
      set_optimizer(mod, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>Tli, "CPX_PARAM_WORKMEM"=>600,"CPX_PARAM_WORKDIR"=> "/home/mariam/julia/julia-1.6.2/bin/CPXsol","CPX_PARAM_NODEFILEIND"=>3))
      Gain=@variable(mod,Gain[l in 1:u])
      x=@variable(mod,x[j in 1:JB,s in 1:M,l in 1:u]>=0,Bin)
      z=@variable(mod,z[b in 1:B,t in 1:T,l in 1:u]>=0,Bin)
    
      q=@variable(mod,q[b in 1:B,t in 1:T,l in 1:u])
      w=@variable(mod,w[b in 1:B,t in 1:T,l in 1:u]>=0,Bin)
      f=@variable(mod,f[j in 1:n,t in 1:T,l in 1:u]>=0)
      Cg=@variable(mod,Cg[ t  in 1:T,l in 1:u]>=0)
      E=@variable(mod,E[b in 1:B,t in 0:T,l in 1:u]>=0)
      y=@variable(mod,y[k in 1:Nbp, t in 0:T, j in 1:JA,l in 1:u]>=0) # température 
      Fin=@variable(mod,Fin[t in 1:T,l in 1:u]>=0)
      Fout=@variable(mod,Fout[t in 1:T,l in 1:u]>=0)
      tot=@variable(mod,tot[t in 1:T,l in 1:u]>=0)
      CC=@variable(mod,CC>=0)
      psi=@variable(mod,psi[j in 1:JA,k in 1:Nbp,t in 1:T,p in 1:pow, l in 1:u], Bin)
    
      #Objective value

      @objective(mod,Max,sum(Gain[l] for l in 1:u))

      @constraint(mod,[ t in 1:T,l in 1:u],sum(Pa[i,j,k,s]*(psi[j,k,t,s,l]-psi[j,k,t,s+1,l]) for j in 1:JA,s in 1:pow-1, k in 1:Nbp if JobA[i,j,k,l]>0)+sum(Pa[i,j,k,pow]*psi[j,k,t,pow,l] for j in 1:JA, k in 1:Nbp)+sum(Pb[i,t,j,s,l]*x[j,s,l] for j in 1:JB,s in 1:M )+sum((q[b,t,l]/dist) for b in 1:B)+p[i,t,l]==L*prod1[i,t,l]+Fin[t,l]+Cg[t,l]-Fout[t,l] )
      @constraint(mod,CC==dist*sum(Cg[t,l] for t in 1:T, l in 1:u))
      
      @constraint(mod,[j in 1:JB,l in 1:u],sum(x[j,s,l] for s in 1:M)==1) #b
      #Constraint with F
      @constraint(mod,[t in 1:T,l in 1:u],Fin[t,l]<=Fbar[t,l])

      @constraint(mod,[t in 1:T,l in 1:u],Fout[t,l]<=P_sous[i])
      @constraint(mod,[t in 1:T,l in 1:u],Fin[t,l]+Cg[t,l]<=P_sous[i])
   
      #Fortz Mars 16th, 2022
      @constraint(mod,[j in 1:JA,k in 1:Nbp, t in 1:T, p in 2:pow, l in 1:u],psi[j,k,t,p,l]<=psi[j,k,t,p-1,l])
      @constraint(mod,[j in 1:JA,k in 1:Nbp, t in 1:T, l in 1:u; JobA[i,j,k,l]>0],psi[j,k,t,1,l]==1)
      
      @constraint(mod,[k in 1:Nbp,t in 1:T,l in 1:u], y[k,t,1,l]==y[k,t-1,1,l]+(delta/Cm[k,i])*(1000*(sum((Pa[i,1,k,p]-Pa[i,1,k,p-1])*psi[1,k,t,p,l] for p in 2:pow )+Pa[i,1,k,1]*psi[1,k,t,1,l])-U[k,i]*(y[k,t-1,1,l]-Yout[l,t]))) 
      @constraint(mod,[k in 1:Nbp], y[k,0,1,1]==y0[k,i])
      @constraint(mod,[k in 1:Nbp,l in 2:u; u>=2],y[k,0,1,l]==y[k,T,1,l-1])
      @constraint(mod,[ k in 1:Nbp,l in 1:u, t in td[k,i,l]:tf[k,i,l]; JobA[i,1,k,l]>0], tmin[k,i,l]<=y[k,t,1,l]<= tmax[k,i,l])

      @constraint(mod,[k in 1:Nbp,l in 1:u,  t in 1:T;t!=tt], y[k,t,2,l]==(y[k,t-1,2,l]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+ ( sum((Pa[i,2,k,p]-Pa[i,2,k,p-1])*psi[2,k,t,p,l] for p in 2:pow )+Pa[i,2,k,1]*psi[2,k,t,1,l])*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) 
      @constraint(mod,[k in 1:Nbp,l in 1:u], y[k,tt,2,l]==ti[k,i,l])
      @constraint(mod,[k in 1:Nbp], y[k,0,2,1]==ti[k,i,1])
      @constraint(mod,[k in 1:Nbp,l in 2:u; u>=2], y[k,0,2,l]==y[k,T,2,l-1])
      @constraint(mod,[ k in 1:Nbp,l in 1:u ; JobA[i,2,k,l]>0], temp_min<=y[k,ec,2,l]<= temp_max)

      @constraint(mod,[b in 1:B,t in 1:T,l in 1:u],dist*d*Pd[i,b]*(z[b,t,l]-1)<=q[b,t,l])
      @constraint(mod,[b in 1:B,t in 1:T,l in 1:u],q[b,t,l]<=dist*c*Pc[i,b]*z[b,t,l])
      @constraint(mod,[b in 1:B,t in 1:T-1,l in 1:u],z[b,t+1,l]-z[b,t,l]<=w[b,t,l])
      @constraint(mod,[b in 1:B,l in 1:u],sum(w[b,t,l] for t in 1:T)<=fi)
      @constraint(mod,[b in 1:B,t in 1:T,l in 1:u],E[b,t,l]==eta*E[b,t-1,l]+q[b,t,l])
      @constraint(mod,[b in 1:B],E[b,0,1]==xi[i,b])
      @constraint(mod,[b in 1:B,l in 2:u;u>=2],E[b,0,l]==E[b,T,l-1])
      @constraint(mod,[t in 1:T,b in 1:B,l in 1:u],E[b,t,l]<=gamma[i,b])

      @constraint(mod,[l in 1:u],Gain[l]==dist*(sum(Fout[t,l]*Pvc[t] for t in 1:T)-sum(Fin[t,l]*Pac[t] for t in 1:T)-Pedf*sum(Cg[t,l] for t in 1:T)))                # Gain après     
    
      # Total consumption en kWh
      @constraint(mod,[t in 1:T,l in 1:u],tot[t,l]==sum(Pa[i,j,k,s]*(psi[j,k,t,s,l]-psi[j,k,t,s+1,l]) for j in 1:JA,s in 1:pow-1, k in 1:Nbp if JobA[i,j,k,l]>0)+sum(Pb[i,t,j,s,l]*x[j,s,l] for j in 1:JB,s in 1:M )+p[i,t,l])
      
     optimize!(mod)   
     return round(JuMP.value(CC),digits=2),round(MOI.get(mod, MOI.SolveTimeSec()), digits=2),round(MOI.get(mod, MOI.RelativeGap()),digits=4),JuMP.value.(Fout),JuMP.value.(Fin),JuMP.value.(tot),JuMP.value.(Gain),JuMP.value.(Cg),JuMP.value.(E),JuMP.value.(q)



end
 
 
 

