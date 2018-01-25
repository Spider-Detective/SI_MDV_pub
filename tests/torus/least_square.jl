function lsf(DR, SVD)
  n_Z = length(DR) ;
  A12 = 0 ; A22 = 0 ; b1 = 0 ; b2 = 0 ;
  for i = 1:n_Z
      dij = DR[i] ; rij = SVD[i] ;
      A12 = A12 + 1./(dij^2)  ;
      A22 = A22 + 1./(dij^4)  ;
      b1  = b1  + rij         ;
      b2  = b2  + rij/(dij^2) ;
  end
  A = [ n_Z A12
	A12 A22 ] ;
  b12 = [ b1
	  b2 ] ;
  abvec = A \ b12
  b = abvec[1] ;
  a = abvec[2] ;

  resvec = zeros(n_Z,1) ;
  zvec = zeros(n_Z,1) ;
  for i = 1:n_Z
      dij = DR[i] ;
      rij = SVD[i] ;
      zij = b + a/(dij^2) ;
      zvec[i,1] = zij ;
      resvec[i,1] = rij - zij ;
  end 
  return zvec, resvec
end
