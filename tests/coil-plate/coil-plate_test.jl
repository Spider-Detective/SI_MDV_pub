using Distances
using PyPlot

include("../../src/funcs.jl")
input = readdlm("d3p_nodeDomainDecomp.orig");

domNum = 64

dRatio = zeros(0)
separation = zeros(0)
farPairs = Int64[]     # Record the pairs with dRatio >=1
n1 = Int64[]           # Size of subdomains
n2 = Int64[]           # Size of subdomains
pts_md = Int64[]      # Number of points chosen in maximally dispersed sets

rank_SI = Int64[]     # Final rank by SI
rank_rrqr = Int64[]   # Final rank by full RRQR
rank_svd = Int64[]    # Final rank by SVD
time_SI = zeros(0)     # Total time of random SI 
time_rrqr = zeros(0)   # Total time of full RRQR 
time_diag = zeros(0)  # Total time of constructing diagonal sub-blocks

seed = 2017
msglvl = 0

tol = 1e-6
kernel = (x,y) -> 1./distance_2(x,y)
useMD = true

i = 1
for ind1 in 1:1
    for ind2 in (ind1+1):domNum 
        println("Calculating pair #", i," (",ind1,", ",ind2,")")

        points1, domain1 = readNodes(input,ind1)         # points1: coordinate of points in domain 1
        points2, domain2 = readNodes(input,ind2)         # domain1: node ID of points in domain 1
        push!(n1, size(points1,2))
        push!(n2, size(points2,2))
        
        msglvl = 0 
        seed = 2017
        
        # Calculate distance ratio
        cent1 = centroid(points1)
        cent2 = centroid(points2)
        maxDist1 = maximum(colwise(Euclidean(), points1, repmat(cent1,1,length(domain1))))
        maxDist2 = maximum(colwise(Euclidean(), points2, repmat(cent2,1,length(domain2))))
        centDist = evaluate(Euclidean(), cent1, cent2)
        distRatio = round(centDist / (maxDist1 + maxDist2),3)
        push!(dRatio, distRatio)
        println("Distance ratio: " ,distRatio)
        
        # Compute the number of MDV
        ncenter = Int(ceil(distRatio^(-2) * 8.48) + 70)
        n_L = size(points2,2)
        ncenter = min(ncenter,n_L)
        
        t1 = time_ns()
        centers, dist, d = findcenters(n_L, points2[1,:], points2[2,:], points2[3,:], ncenter, seed, msglvl); 
        push!(pts_md, ncenter)
        
        # Do SI
        XYb = meshKernelFull(kernel, points1, points2[:,centers]);
        Q1, R1, P1, r1 = rrqr(XYb',tol)
        Xh = points1[:,P1[1:r1]]
        XhY = meshKernelFull(kernel, Xh, points2);
        Q2, R2, P2, r2 = rrqr(XhY,tol)
        nh = max(r1,r2)
        Xh = points1[:,P1[1:nh]]
        Yh = points2[:,P2[1:nh]]
        U = meshKernelFull(kernel, points1, Yh);
        # In this case, the factorization is (U * Sinv * V') = U * (S \ V')
        S = lufact(meshKernelFull(kernel, Xh, Yh)); 
        V = meshKernelFull(kernel, Xh, points2)';
        t2 = time_ns()
        println("Total time of SI: ", (t2-t1)/1e9)
        println("Final rank of SI: ", nh)
        push!(time_SI, (t2-t1)/1e9)
        push!(rank_SI, nh)
        
        # Do full RRQR
        t1 = time_ns()
        Atrue = meshKernelFull(kernel, points1, points2)
        Q, R, P, r = rrqr(Atrue, tol)
        t2 = time_ns()
        println("Total time of full RRQR: ", (t2-t1)/1e9)
        #println("Final rank of full RRQR: ", r)
        push!(time_rrqr, (t2-t1)/1e9)
        push!(rank_rrqr, r)
        
        # Do SVD
        t1 = time_ns()
        Atrue = meshKernelFull(kernel, points1, points2)
        Acopy=copy(Atrue)
        r_svd = rank_eps_fro(Acopy,tol)[1]
        t2 = time_ns()
        println("Final rank of SVD: ", r_svd)
        push!(rank_svd, r_svd)

        if distRatio>=1
            push!(farPairs, i)
        end
        println()
        i = i + 1
    end
end
# Construct all diagonal sub-blocks
i = 1
for ind in 1:domNum
        println("Calculating pair #", i," (",ind,", ",ind,")")
        points1, domain1 = readNodes(input,ind)         # points1: coordinate of points in domain 1
       
        t1 = time_ns()
        Atrue = meshKernelFull(kernel, points1, points1)
        t2 = time_ns()
        println("Total time of construction: ", (t2-t1)/1e9)
        push!(time_diag, (t2-t1)/1e9)
        println()
        i = i + 1
end

writeIn = false
if writeIn
  # Save all data for each case
  open("MD_data_full_1_new.txt", "w") do f
      write(f, "Original decomposition of EM coil. Total subdomain N = 64. \n")
      write(f, "All pairs computed. tol = 1e-6. Use Frobenius norm to cutoff. \n") 
      write(f, "Use maximally dispersed set to pick up Ybar and do SI.\n")
      write(f, "distRatio = ||Centroid_{i,j}|| / (Radius_i + Radius_j) \n")
      write(f, "Choose the points in dispersed set as ncenter = Int(ceil(distRatio^(-2) * 8.48 + 70)).\n")
      write(f, "Total # of pairs: ")
      writedlm(f, length(dRatio))
      write(f, "Distance ratio: \n")
      writedlm(f, dRatio')
      write(f, "Size of subdomain n1: \n")
      writedlm(f, n1')
      write(f, "Size of subdomain n2: \n")
      writedlm(f, n2')
      write(f, "# of points in maximally dispersed sets: \n")
      writedlm(f, pts_md')
      write(f, "Rank by SI: \n")
      writedlm(f, rank_SI')
      write(f, "Rank by full RRQR: \n")
      writedlm(f, rank_rrqr')
      write(f, "Rank by SVD: \n")
      writedlm(f, rank_svd')
      write(f, "Time of SI: \n")
      writedlm(f, time_SI')
      write(f, "Time of full RRQR: \n")
      writedlm(f, time_rrqr')
      write(f, "Time for constructing diagonal sub-blocks: \n")
      writedlm(f, time_diag')
  end
end
