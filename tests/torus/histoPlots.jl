using PyPlot

data = readdlm("MD_data_full_1.txt")
rSI1 = data[15,1:2016]
rRRQR1 = data[17,1:2016]
rSVD1 = data[19,1:2016]
data = readdlm("MD_data_full_2.txt")
rSI2 = data[15,1:2016]
rRRQR2 = data[17,1:2016]
rSVD2 = data[19,1:2016]
data = readdlm("MD_data_full_3.txt")
rSI3 = data[15,1:2016]
rRRQR3 = data[17,1:2016]
DR = data[8,1:2016]
rSVD3 = data[19,1:2016];

fp = find(DR->DR.>1,DR)
rank_diff1 = rSI1 .- rRRQR1
rank_diff2 = rSI2 .- rRRQR2
rank_diff3 = rSI3 .- rRRQR3

x = rank_diff1[fp]
hist1 = Int64[]
range1 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range1
    push!(hist1, length(find(x->x==i,x)))
end

x = rank_diff2[fp]
hist2 = Int64[]
range2 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range2
    push!(hist2, length(find(x->x==i,x)))
end

x = rank_diff3[fp]
hist3 = Int64[]
range3 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range3
    push!(hist3, length(find(x->x==i,x)))
end
figure(1)
plot(range1,hist1,"ob--""",label="a=70,b=8.48")
plot(range2,hist2,"oy--",label="a=80,b=8.48")
plot(range3,hist3,"og--",label="a=90,b=8.48")
legend(bbox_to_anchor=[0.3,1],loc=0,borderaxespad=0)
xlabel("\$r_{SI-MDV}-r_{rrqr}, DR>1\$")
ylabel("Counts")
grid("on")
savefig("torus_rrqr_SI-MDV_far.pdf")

rank_diff1 = rSI1 .- rSVD1
rank_diff2 = rSI2 .- rSVD2
rank_diff3 = rSI3 .- rSVD3
x = rank_diff1[fp]
hist1 = Int64[]
range1 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range1
    push!(hist1, length(find(x->x==i,x)))
end

x = rank_diff2[fp]
hist2 = Int64[]
range2 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range2
    push!(hist2, length(find(x->x==i,x)))
end

x = rank_diff3[fp]
hist3 = Int64[]
range3 = collect(linspace(
                        minimum(x),
                        maximum(x),
                        maximum(x)-minimum(x)+1
                        ))
for i in range3
    push!(hist3, length(find(x->x==i,x)))
end

figure(2)
plot(range1,hist1,"ob--""",label="a=70,b=8.48")
plot(range2,hist2,"oy--",label="a=80,b=8.48")
plot(range3,hist3,"og--",label="a=90,b=8.48")
legend(bbox_to_anchor=[0.3,1],loc=0,borderaxespad=0)
xlabel("\$r_{SI-MDV}-r_{SVD}, DR>1\$")
ylabel("Counts")
grid("on")
savefig("torus_svd_SI-MDV_far.pdf")
