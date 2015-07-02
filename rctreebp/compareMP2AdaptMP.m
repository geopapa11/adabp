E = [1  2;
     1  3;
     5  6;
     3  4;
     3  5];
  
phi{1} = [.8 .2]';  phi{2} = [.4 .6]';  phi{3} = [.6 .1 .3]';  phi{4} = [.1 .9]';
phi{5} = [.4 .6]';  phi{6} = [.3 .1 .2 .4 .0]';

psi{1} = [.9 .1; 
          .3 .7];

psi{2} = [.1 .7 .2;
          .2 .3 .5];

psi{3} = [.3 .2 .1 .2 .1; 
          .1 .6 .1 .1 .1];

psi{4} = [.2 .8;
          .7 .3;
          .6 .4];

psi{5} = [.7 .3;
          .4 .6;
          .2 .8];
      
N  = length(phi);
Ne = size(E,1);

mapseq1 = bp_maxsum(E, phi, psi);

% Build Adjacent (sparse) matrix - rows are variables, cols are factors
% First N factors correpond to node potentials, remaining Ne to pairwise
% potentials
A = sparse([1:N, E(:)'], [1:N, repmat(N+1:N+Ne,1,2)], 1);

% Build factor functions
F = [phi, psi];

% Build RC-tree structure for Max-Product
T = rctreeMP(A,F,8);
mapseq2 = double(T.map);

ibp = IncrBP_NEW(E,phi,psi,'max');
mapseq3 = ibp.mapseq();

mapseq1-mapseq2
mapseq1-mapseq3

L1 = loglik(mapseq1, E, phi, psi)
L2 = loglik(mapseq2, E, phi, psi)
L3 = loglik(mapseq3, E, phi, psi)

disp('Incorporating a new measurement at node 5');
disp('-----------------------------------------');
j = 5;
phi{j} = [.1 .9]';
mapseq1 = bp_maxsum(E, phi, psi);
T       = updateFactorMP(T,j,phi{j});
mapseq2 = double(T.map);

ibp.setNodePot(j, phi{j});
mapseq3 = ibp.mapseq();

mapseq1-mapseq2
mapseq1-mapseq3

L1 = loglik(mapseq1, E, phi, psi)
L2 = loglik(mapseq2, E, phi, psi)
L3 = loglik(mapseq3, E, phi, psi)


disp('Incorporating a new measurement at node 1');
disp('-----------------------------------------');
j_prv = j;
j = 1;
phi{j} = [.1 .9]';
mapseq1 = bp_maxsum(E, phi, psi);
T       = updateFactorMP(T,j,phi{j});
mapseq2 = double(T.map);

ibp.setNodePot(j, phi{j});
ibp.propagate(j_prv, j, true);
mapseq3 = ibp.mapseq();

mapseq1-mapseq2
mapseq1-mapseq3

L1 = loglik(mapseq1, E, phi, psi)
L2 = loglik(mapseq2, E, phi, psi)
L3 = loglik(mapseq3, E, phi, psi)


disp('Incorporating a new measurement at node 3');
disp('-----------------------------------------');
j_prv = j;
j = 3;
phi{j} = [.2 .5 .3]';
mapseq1 = bp_maxsum(E, phi, psi);
T       = updateFactorMP(T,j,phi{j});
mapseq2 = double(T.map);

ibp.setNodePot(j, phi{j});
ibp.propagate(j_prv, j, true);
mapseq3 = ibp.mapseq();

mapseq1-mapseq2
mapseq1-mapseq3

L1 = loglik(mapseq1, E, phi, psi)
L2 = loglik(mapseq2, E, phi, psi)
L3 = loglik(mapseq3, E, phi, psi)
