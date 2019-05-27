function [ A ] = gen_mat( problem )
%gen_mat(problem) : generates a symmetric positive definite matrix of size n and
%norm equal to one, either according to the Wishart distribution or with a
%specific condition number

% problem is a structure containing the fields : 
%        problem.size
%        (problem.cond)
%        (problem.version)

% Author: E.Massart

if ~isfield(problem,'cond')
    A = randn(problem.size,problem.size);
    A = A'*A;
else
    if ~isfield(problem, 'version')
        problem.version = 1;
    end
    
    if problem.version == 1
        %     disp('version1');
        [Q,~] = qr(rand(problem.size));
        D = diag([rand(1,problem.size-1)+1,10^(-problem.cond)]);
        A = Q*D*Q';
    elseif problem.version == 2
        %     disp('version2');
        [Q,~] = qr(rand(problem.size));
        D = diag([1,10^(-problem.cond)+10^(-problem.cond)*rand(1,problem.size-1)]);
        A = Q*D*Q';
    elseif problem.version == 3
        %     disp('version3');
        [Q,~] = qr(rand(problem.size));
        if mod(problem.size,2) == 0
            D = diag([1+rand(1,problem.size/2),10^(-problem.cond)+10^(-problem.cond)*rand(1,problem.size/2)]);
        else
            D = diag([1+rand(1,floor(problem.size/2)+1),10^(-problem.cond)+10^(-problem.cond)*rand(1,floor(problem.size/2))]);
        end
        A = Q*D*Q';
    end
end

nor = norm(A,'fro');
A = A./nor;

lambdas = eig(A);
e = find(lambdas<=0);

if(size(e)~= 0)
    disp('The matrix is not positive definite');
end

end




