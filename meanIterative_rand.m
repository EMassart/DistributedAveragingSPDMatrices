function [ M_fin, info] = meanIterative_rand( A, options )
%MEANITERATIVE_RAND(A) estimates the Karcher mean using an algorithm based on
%the Iterative algorithm proposed in [1].
%The cycle along which the matrices are averaged is randomly chosen at each
%iteration.
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

    G = randperm(N);

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





