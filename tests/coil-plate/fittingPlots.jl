using PyPlot
include("least_square.jl")

data = readdlm("MD_data_4.txt")
rSI3 = data[18,1:2016]
DR3 = data[8,1:2016]
rSVD3 = data[22,1:2016];

figure(1)
plot(DR3,rSVD3,".b",label="SVD")
plot(DR3,rSI3,".r",label="SI-MDV,a=90,b=8.44")
legend(bbox_to_anchor=[1,1],loc=0,borderaxespad=0)
grid("on")
savefig("torus_SVD_SI-MDV.pdf")

figure(2)
zvec, resvec, a, b = lsf(DR3, rSVD3)
plot(DR3, rSVD3, ".b",label="SVD") ;
sortD = sortperm(DR3)
plot(DR3[sortD], zvec[sortD], "r",label="Least Square Fitting") ;
legend(bbox_to_anchor=[1,1],loc=0,borderaxespad=0)
grid("on") ;
xlabel("distance ratio") ;
ylabel("rank") ;
println("a=",round(a,2))
println("b=",round(b,2))
savefig("torus_SVD_fitting.pdf")

figure(3)
fp3 = find(DR3->DR3.>1,DR3)
DR3far = DR3[fp3]
rSVD3far = rSVD3[fp3]
zvec, resvec, a, b = lsf(DR3far, rSVD3far)
plot(DR3far, rSVD3far, ".b",label="SVD") ;
sortD = sortperm(DR3far)
plot(DR3far[sortD], zvec[sortD], "r",label="Least Square Fitting") ;
legend(bbox_to_anchor=[1,1],loc=0,borderaxespad=0)
grid("on") ;
xlabel("distance ratio") ;
ylabel("rank") ;
println("a=",round(a,2))
println("b=",round(b,2))
savefig("torus_SVD_fitting_far.pdf")
