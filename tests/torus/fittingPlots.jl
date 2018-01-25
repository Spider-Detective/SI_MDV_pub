using PyPlot
include("least_square.jl")

data = readdlm("MD_data_full_1.txt")
rSI1 = data[15,1:2016]
rRRQR1 = data[17,1:2016]
rSVD1 = data[19,1:2016]
n1 = data[10,1:2016]
n2 = data[12,1:2016]
DR1 = data[8,1:2016]
md1 = data[13,1:2016]
tSI1 = data[21,1:2016]
tRRQR1 = data[23,1:2016];

data = readdlm("MD_data_full_2.txt")
rSI2 = data[15,1:2016]
rRRQR2 = data[17,1:2016]
DR2 = data[8,1:2016]
md2 = data[13,1:2016]
rSVD2 = data[19,1:2016]
tSI2 = data[21,1:2016]
tRRQR2 = data[23,1:2016];

data = readdlm("MD_data_full_3.txt")
rSI3 = data[15,1:2016]
rRRQR3 = data[17,1:2016]
DR3 = data[8,1:2016]
md3 = data[13,1:2016]
rSVD3 = data[19,1:2016];
tSI3 = data[21,1:2016]
tRRQR3 = data[23,1:2016];

figure(1)
plot(DR3,rSVD3,".b",label="SVD")
plot(DR3,rSI3,".r",label="SI-MDV,a=90,b=8.44")
legend(bbox_to_anchor=[1,1],loc=0,borderaxespad=0)
grid("on")
savefig("torus_SVD_SI-MDV.pdf")

figure(2)
zvec, resvec = lsf(DR3, rSVD3)
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
zvec, resvec = lsf(DR3far, rSVD3far)
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
