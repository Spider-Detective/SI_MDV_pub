function distance_2{T}(x::Array{T, 1}, y::Array{T, 1})
    @assert size(x) == size(y)
    return sqrt(sum((x - y).^2,1))
end

function distance_2{T}(x::Array{T, 2}, y::Array{T, 2})
    @assert size(x) == size(y)
    return sqrt(sum((x - y).^2,1))[:]
end

# Builds the matrix entry by entry
function meshKernelFull(kernel, evalPts1::Array{Float64, 2}, evalPts2::Array{Float64, 2})
    n1 = size(evalPts1)[2]
    n2 = size(evalPts2)[2]

    X = zeros(Float64, size(evalPts1)[1], n1*n2)
    Y = zeros(Float64, size(evalPts2)[1], n1*n2)
    k = 1
    for j = 1:n2 # Careful with the ordering
        for i = 1:n1
            X[:,k] = evalPts1[:,i]
            Y[:,k] = evalPts2[:,j]
            k = k+1;
        end
    end
    A = reshape(kernel(X, Y), (n1, n2))

    return A
end
# A method for calculating the centroid of a bunch of points
function centroid(pts)
    dim = size(pts,1)
    npts = size(pts,2)
    cent = zeros(dim, 1)
    for i in 1:dim
        cent[i] = sum(pts[i,:]) / npts
    end
    return cent
end

# Read in all the points in a pair of subdomains
function readNodes(input, index)
    domain = zeros(0)
    # Read the subdomain nodes
    for i in 1:length(input)
        if input[i] == "*SET_NODE"
            if input[i+1] == index
                j = i+2     # skip 2 lines of headings
                while input[j] != "*SET_NODE" && j <= size(input,1)
                    push!(domain, input[j])
                    j += 1
                end
            end
        end
    end
    # Get the coordinates of the subdomain nodes, for 3D
    points = zeros(length(domain),3)
    for i in 1:length(domain)
        pts = input[Int(domain[i])+2,2:4]
        for j in 1:2
            s = pts[j]
            pts[j] = parse(Float64,s[1:(end-1)])
            points[i,j] = pts[j]
            end
            points[i,3] = pts[3]
        end
        points = points'
    return points, domain
end

function cut_spectrum_fro(s, tol; bypasssort=false)
    if ! bypasssort
        @assert issorted(s, rev=true)
    else
        if ! issorted(s, rev=true)
            warn(string("Not sorted array : ", string(s)))
        end
    end
    @assert all(s .> 0)
    n  = length(s)
    ns = vecnorm(s)
    for j = 0:n-1
        if vecnorm(s[j+1:end]) <= tol #* ns
            return j
        end
    end  
    return n
end

# perform SVD
function rank_eps_fro(A, tols)
    (U,s,V) = svd(A)
    @assert issorted(s, rev=true)
    ranks = Array{Int64, 1}(length(tols))
    for i = 1:length(tols)
        j = cut_spectrum_fro(s, tols[i])
        ranks[i] = j
        # Actually check accuracy
        #err = vecnorm(U[:,1:j]*diagm(s[1:j])*V[:,1:j]' - A)/vecnorm(A)
    end
    return ranks
end
# column pivoting RRQR
function rrqr(A,tol)
    m, n = size(A)
    Q = zeros(m,n)
    R = zeros(m,n)
    rank = 0
    ind = collect(1:n)
    Acopy = copy(A)
    
    nrm = vecnorm(A)
    for r in 1:min(m,n)
        cn = zeros(n,1)
        for j=r:n
            cn[j] = vecnorm(A[:,j])
        end
        tau = maximum(cn)
        if tau!=0
            k = indmax(cn)
        end
        if k>r
            tmp = ind[r]; ind[r] = ind[k]; ind[k] = tmp;
            tmp = A[:,r]; A[:,r] = A[:,k]; A[:,k] = tmp;
            tmp = R[1:r-1,r]; R[1:r-1,r] = R[1:r-1,k]; R[1:r-1,k] = tmp;
        end
        R[r,r] = tau
        
        #if (R[1,1]<=tol || R[r,r]/R[1,1]<=tol)
        if (R[1,1]<=tol || vecnorm(triu(A[r+1:end,r+1:end]))<=tol)
        #if (R[1,1]<=tol || vecnorm(triu(A[r+1:end,r+1:end]))/nrm<=tol) # Golub van Loan section 5.5.7 eq 5.5.6
            r -= 1
            rank = r
            break;
        end
        Q[:,r] = A[:,r]/R[r,r];
        R[r,r+1:n] = Q[:,r]' * A[:,r+1:n]
        A[:,r+1:n] = A[:,r+1:n] - Q[:,r]*R[r,r+1:n]'
        rank = r
    end
    
    if rank==0
        Q = []
        R = []
    else
        Q = Q[:,1:rank]
        R = R[1:rank,:]
    end
    return Q, R, ind, rank
end

# function [ centers, dist, d ] ...
#                 = findcenters (npts, x, y, z, ncenter, seed, msglvl )
# input --
#     npts -- # of points
#     x[]  -- x locations, size npts
#     y[]  -- y locations, size npts
#     z[]  -- z locations, size npts
#     ncenter -- # of desired centers, maximally dispersed vertices
#     seed    -- random number seed
#     msglvl  -- message level
# output --
#     centers[] -- vector of vertex ids, size ncenter
#     dist[]    -- matrix of distances, size npts x ncenter
#     d[]       -- vector of minimum distances to previous centers, 
#                  size ncenter
# created -- a long time ago, cca
function findcenters(npts, x, y, z, ncenter, seed, msglvl)
    srand(seed) 
    dist = zeros(npts, ncenter) 
    # pick a node at random
    jj = Int(ceil(npts*rand())) 
    for ii = 1:npts
       dist[ii,1] = sqrt( (x[ii] - x[jj])^2 + (y[ii] - y[jj])^2 + (z[ii] - z[jj])^2 ) 
    end
    # find node "center" that is the maximum distance from jj
    center = 0 
    maxdist = 0 
    for ii = 1:npts
       if maxdist < dist[ii,1]
          maxdist = dist[ii,1] 
          center = ii 
       end
    end
    #  use minimum distance to a previous center
    centers = zeros(Int64, ncenter) 
    centers[1] = center 
    for jj = 1:ncenter
       if msglvl > 2
          fprintf("\n jj = %d", jj) 
       end
       center = centers[jj]
       if msglvl > 2
          fprintf(", center = %d", center) 
       end
       for ii = 1:npts
          dist[ii,jj] = sqrt( (x[ii] - x[center])^2  + (y[ii] - y[center])^2 + (z[ii] - z[center])^2 )
       end
       maxdist = 0 
       next = 0
       for ii = 1:npts
          mindist = dist[ii,1] 
          for kk = 2:jj
             if mindist > dist[ii,kk]
                mindist = dist[ii,kk] 
             end
          end
          if maxdist == mindist
             if next == 0
                next = ii 
             end
          elseif maxdist < mindist
             maxdist = mindist 
             next = ii 
          end
       end
       if jj == ncenter
          break 
       end
       centers[jj+1,1] = next 
    end
    
    # compute the minimal distances d(k) = min(d(i,j)) for (i,j) in C_k
    d = zeros(ncenter,1) 
    d[1] = dist[1,1] 
    if msglvl > 2
       fprintf("\n d(%d) = %.6f", 1, d[1]) 
    end
    for k = 2:ncenter
       K = centers[k] 
       if msglvl > 2
          fprintf("\n k %d, K %d", k, K) 
       end
       J = centers[1]
       mindist = dist[J, k] 
       if msglvl > 2
          fprintf("\n    dist(%d,%d) = %.6f", J, K, dist[J,k]) 
       end
       for jj = 2:k-1
          J = centers[jj]
          if msglvl > 2
             fprintf("\n    dist(%d,%d) = %.6f", J, K, dist[J,k]) 
          end
          mindist = min(mindist, dist[J,k]) 
       end
       d[k] = mindist 
       if msglvl > 2
          fprintf("\n d(%d) = %.6f", k, d[k]) 
       end
    end
    d[1] = 0 

    return centers, dist, d
end

