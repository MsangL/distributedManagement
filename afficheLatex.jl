function affiche(inj,stock,rest,need,prod1,conso,qP,share,ro)
         s=""
         s1=""
         s2=""
         s3=""
         s4=""
         s5=""
         s6=""
         s7=""
         for t in 1:T
                  s=s*"("*string(t)*","*string(inj[t])*")"
                  s1=s1*"("*string(t)*","*string(stock[t])*")"
                  s2=s2*"("*string(t)*","*string(rest[t])*")"
                  s3=s3*"("*string(t)*","*string(need[t])*")"
                  s4=s4*"("*string(t)*","*string(prod1[t])*")"
                  s5=s5*"("*string(t)*","*string(conso[t])*")"
                  s6=s6*"("*string(t)*","*string(qP[t])*")"
                  s7=s7*"("*string(t)*","*string(share[t])*")"
         end
    open("Coube","a") do io
          println(io,"
              \\begin{subfigure}[t]{.45\\linewidth}
               \\centering
               \\begin{tikzpicture}[scale=0.8]
               \\begin{axis}[grid= major ,
                xlabel = {Time} ,
                ylabel = {} ,
                xmin = 0, xmax = 24,
                ymin = 0, ymax = 18,
                legend entries={Total cons,Total need,Injected,Shared,Charge BSS,rest},
                legend style={at={(0,1)},anchor=north west}]")
          #println(io,"\\addplot[green, very thick,mark=*] coordinates{",s4,"};")
          println(io,"\\addplot[gray, very thick,mark=o] coordinates{",s5,"};")
          println(io,"\\addplot[orange, very thick,mark=+] coordinates{",s3,"};")
          println(io,"\\addplot[red, very thick,mark=x] coordinates{",s,"};")
          #println(io,"\\addplot[brown, very thick,mark=|] coordinates{",s1,"};")
          println(io,"\\addplot[green, very thick,mark=o] coordinates{",s7,"};")
          println(io,"\\addplot[blue, very thick,mark=o] coordinates{",s6,"};")
          println(io,"\\addplot[blue, very thick,mark=o] coordinates{",s2,"};")
          println(io,"  \\end{axis}
          \\end{tikzpicture}    
          \\caption{ Distributed \$\\rho\$=",ro,".}
            \\end{subfigure}
            \\hfill")
     end
end

               # legend entries={Production(kW),Conso,Besoin totale(kWh),Injectée(kWh),Stockage(kWh),Charge BSS, Reste(kWh)},
