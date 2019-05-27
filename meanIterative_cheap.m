function [ M_fin, info] = meanIterative_cheap( A, options )
%MEANITERATIVE_CHEAP(A) estimates the Karcher mean using an algorithm based on
%the Iterative algorithm proposed in [1].
%The cycle is chosen at each iteration in order to maximize the
%communication speed between the nodes.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

%References : 
%[1] Miklós Pálfia. A multivariable extension of two-variable matrix means. SIAM Journal
%on Matrix Analysis and Applications, 32(2):385–393, 2011.

tstart = cputime;
tcumulstop = 0;   
    
N = length(A);
C = zeros(N,N);                                     %will contain the weights of the influence of each input matrix in the current nodes
                                                    %(used to compute the new cycle at each iteration)
nIterMax = options.n_iter;
time = zeros(1,nIterMax+1);
d_min = zeros(1,nIterMax+1);
d_max = zeros(1,nIterMax+1);
d_mean = zeros(1,nIterMax+1);
iter = 0;
count = 0;                                          %records the number of means of two matrices evaluated
M = cell(1,nIterMax);

tcumulstart = cputime;
d = zeros(1,N);
for i = 1:N
    d(i) = geo_dist(A{i},options.ref);
end
d_ref = mean(d);
d_min(1) = min(d);
d_max(1) = max(d);
d_mean(1) = mean(d);
tcumulstop = tcumulstop + cputime - tcumulstart;

while(iter<nIterMax)
     
    iter = iter+1;
    %computation of the cycle (the assumption is made here that the
    %sequence of cycles is precomputed, therefore the required time is not recorded)
    if (iter==1)
        
        G = 1:N;                   
        for i = 1:N-1                       %C is firstly defined as the transposed of the incidence matrix
            C(i,G(i)) = 1;
            C(i,G(i+1)) = 1;
        end
        C(N,G(N)) = 1;
        C(N,G(1)) = 1;
        
    else
        
        D_cycle = triu(squareform(pdist(C)));
        G = idealMap(D_cycle);                      %computes the new cycle
        Cbis = C;                                   %updates the matrix C 
        for i = 1:N-1
            C(i,:) = (Cbis(G(i),:) + Cbis(G(i+1),:))./2;  
        end
        C(N,:) = (Cbis(G(N),:) + Cbis(G(1),:))./2;
    end
    
    B = A;
    for i = 1:N-1
        A{i} = mean2(B{G(i)}, B{G(i+1)});
    end
    A{N} = mean2(B{G(N)}, B{G(1)});
    time(iter+1) = cputime - tstart - tcumulstop;
    
    tcumulstart = cputime;
    d = zeros(1,N);
    for i = 1:N
        d(i) = geo_dist(A{i},options.ref);
    end
    d_min(iter+1) = min(d);
    d_max(iter+1) = max(d);
    d_mean(iter+1) = mean(d);

    M{iter} = sum(cat(3,A{:}),3)./N;     
    tcumulstop = tcumulstop + cputime - tcumulstart;
        
    count = count + N;

end


M_fin = M{end};
info.M_rec = M;
info.d_min = d_min./d_ref;
info.d_max = d_max./d_ref;
info.d_mean = d_mean./d_ref;
info.time = time;
info.t_fin = cputime-tstart;

end





