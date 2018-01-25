data = readdlm("MD_data_2.txt")
rSI1 = data[18,1:2016]
rRRQR1 = data[20,1:2016]
rSVD = data[22,1:2016]
n1 = data[13,1:2016]
n2 = data[15,1:2016]
tSI1 = data[24,1:2016]
tRRQR1 = data[26,1:2016];
time_diag = data[28,1:64]

data = readdlm("MD_data_3.txt")
rSI2 = data[18,1:2016]
rRRQR2 = data[20,1:2016]
rSVD2 = data[22,1:2016]
tSI2 = data[24,1:2016]
tRRQR2 = data[26,1:2016];

data = readdlm("MD_data_4.txt")
rSI3 = data[18,1:2016]
rRRQR3 = data[20,1:2016]
DR = data[8,1:2016]
rSVD3 = data[22,1:2016];
tSI3 = data[24,1:2016]
tRRQR3 = data[26,1:2016];

N = 20794
fp = find(DR->DR.>1,DR) 
np = find(DR->DR.<=1,DR)

# Storage results, all pairs
storage_SI = rSI1[fp].*(n1[fp] .+ n2[fp])
storage_rrqr = rRRQR3[np].*(n1[np] .+ n2[np])
storage_diag = n1[1]^2 + sum(n2[1:63].^2)
stPer = (sum(storage_SI) .+ sum(storage_rrqr) +storage_diag) / N^2
println("Storage percentage (a=70): ", round(stPer,5))
storage_SI = rSI2[fp].*(n1[fp] .+ n2[fp])
storage_rrqr = rRRQR3[np].*(n1[np] .+ n2[np])
stPer = (sum(storage_SI) .+ sum(storage_rrqr) +storage_diag) / N^2
println("Storage percentage (a=80): ", round(stPer,5))
storage_SI = rSI3[fp].*(n1[fp] .+ n2[fp])
storage_rrqr = rRRQR3[np].*(n1[np] .+ n2[np])
stPer = (sum(storage_SI) .+ sum(storage_rrqr) +storage_diag) / N^2
println("Storage percentage (a=90): ", round(stPer,5))
storage_SI = rSVD[fp].*(n1[fp] .+ n2[fp])
storage_rrqr = rRRQR3[np].*(n1[np] .+ n2[np])
stPer = (sum(storage_SI) .+ sum(storage_rrqr) +storage_diag) / N^2
println("Storage percentage (SVD): ", round(stPer,5))
println()
# Timing results
println("Timing results for torus case (DR>1): ")
println("a=70, Time for SI-MDV: ", round(sum(tSI1[fp]),2))
println("a=70, Time for RRQR: ", round(sum(tRRQR3[fp]),2))
println("Percentage: ",round(sum(tSI1[fp]),2)/round(sum(tRRQR3[fp]),2))
println()
println("a=80, Time for SI-MDV: ", round(sum(tSI2[fp]),2))
println("a=80, Time for RRQR: ", round(sum(tRRQR3[fp]),2))
println("Percentage: ",round(sum(tSI2[fp]),2)/round(sum(tRRQR3[fp]),2))
println()
println("a=90, Time for SI-MDV: ", round(sum(tSI3[fp]),2))
println("a=90, Time for RRQR: ", round(sum(tRRQR3[fp]),2))
println("Percentage: ",round(sum(tSI3[fp]),2)/round(sum(tRRQR3[fp]),2))
println()
println("Timing results for torus case (all pairs): ")
println("a=70, Time with SI-MDV: ", round(sum(tSI1[fp])+sum(tRRQR1[np])+sum(time_diag),2))
println("a=70, Time without SI-MDV: ", round(sum(tRRQR3)+sum(time_diag),2))
println("Percentage: ",round(sum(tSI1[fp])+sum(tRRQR3[np])+sum(time_diag),2)/round(sum(tRRQR3)+sum(time_diag),2))
println()
println("a=80, Time with SI-MDV: ", round(sum(tSI2[fp])+sum(tRRQR2[np])+sum(time_diag),2))
println("a=80, Time without SI-MDV: ", round(sum(tRRQR3)+sum(time_diag),2))
println("Percentage: ",round(sum(tSI2[fp])+sum(tRRQR3[np])+sum(time_diag),2)/round(sum(tRRQR3)+sum(time_diag),2))
println()
println("a=90, Time with SI-MDV: ", round(sum(tSI3[fp])+sum(tRRQR3[np])+sum(time_diag),2))
println("a=90, Time without SI-MDV: ", round(sum(tRRQR3)+sum(time_diag),2))
println("Percentage: ",round(sum(tSI3[fp])+sum(tRRQR3[np])+sum(time_diag),2)/round(sum(tRRQR3)+sum(time_diag),2))
