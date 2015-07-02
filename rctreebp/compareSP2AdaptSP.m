E = [1  2;
     1  3;
     5  6;
     3  4];
  
phi{1} = [.8 .2]';  phi{2} = [.4 .6]';  phi{3} = [.6 .1 .3]';  phi{4} = [.1 .9]';
phi{5} = [.4 .6]';  phi{6} = [.3 .1 .2 .1 .2]';

psi{1} = [.9 .1; 
          .3 .7];
psi{2} = [.1 .7 .2; 
          .2 .3 .5];
psi{3} = [.3 .2 .1 .2 .1; 
          .1 .6 .1 .1 .1];
psi{4} = [.2 .8;
          .7 .3;
          .6 .4];

N  = length(phi);
Ne = size(E,1);

b1 = bp(E, phi, psi);


% Build Adjacent (sparse) matrix - rows are variables, cols are factors
% First N factors correpond to node potentials, remaining Ne to pairwise
% potentials
A = sparse([1:N, E(:)'], [1:N, repmat(N+1:N+Ne,1,2)], 1);

% Build factor functions
F = [phi, psi];

% Build RC-tree structure for Sum-Product
T  = rctreeSP(A,F);

b2 = cell(N,1); 
for i = 1:N
    b2{i} = marginal(T,i);
end

for i = 1:N
    foo = b1{i} - b2{i}
end

disp('Incorporating a new measurement at node 5');
disp('-----------------------------------------');
j = 5;
phi{j} = [.9 .1]';
b1 = bp(E, phi, psi);
T  = updateFactorSP(T,j,phi{j});
for i = 1:N
    b2{i} = marginal(T,i);
end

for i = 1:N
    foo = b1{i} - b2{i}
end


disp('Incorporating a new measurement at node 1');
disp('-----------------------------------------');
j = 1;
phi{j} = [.1 .9]';
b1 = bp(E, phi, psi);
T  = updateFactorSP(T,j,phi{j});
for i = 1:N
    b2{i} = marginal(T,i);
end

for i = 1:N
    foo = b1{i} - b2{i}
end


disp('Incorporating a new measurement at node 3');
disp('-----------------------------------------');
j = 3;
phi{j} = [.2 .5 .3]';
b1 = bp(E, phi, psi);
T  = updateFactorSP(T,j,phi{j});
for i = 1:N
    b2{i} = marginal(T,i);
end

for i = 1:N
    foo = b1{i} - b2{i}
end
