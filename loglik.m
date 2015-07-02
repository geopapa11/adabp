function L = loglik(seq, E, phi, psi)
% function L = loglik(seq, E, phi, psi)
% 
% INPUT
%     N:   1x1   scalar         -  Number of nodes (variables)
%    Ne:   1x1   scalar         -  Number of edges
%  card:   1x1   scalar         -  Cardinality (alphabet) of the variables
%                                  That is, each variable x_i takes values 
%                                  from {1, ..., card}, that is, 
%                                  X_i \in \mathcal{X}_i = {1, ..., card_i}
%   seq:  1xN double array      -  Sequence
%                                  Sequence of variable values
%     E: Nex2 double matrix     -  Edges
%                                  Each row corresponds to an edge
%   phi:  Nx1 cell array        -  Node potentials
%         or                       Each element is a card(i) x 1 vector
%         card x 1 double array    Node potential: common across all nodes
%   psi:  Nex1 cell array       -  Pairwise potentials
%                                  Each element represents a pairwise potential 
%                                  between nodes defined by double matrix E
%                                  Each element (i,j) is a card(i) x card(j) vector
%   or  card x card             -  Pairwise potential 
%       double matrix              (common across all edges)
% 
% 
% OUTPUT
%     L:  1x1 scalar            -  Log-likelihood of the sequence
% 
% 
% EXAMPLES
%
% E.g. assume the simple graphical model:
%
%       1
%     /   \
%    2     3
%           \
%            4
% 
% 
% We have the following edges, node and pairwise potentials, respectively:
% 
% seq = [2 3 1 4];
% 
% E = [1  2;
%      1  3; 
%      3  4];
% 
% phi{1} = [.8 .2]';  phi{2} = [.2 .3 .4]';  phi{3} = [.9 .1]';  phi{4} = [.2 .1 .6 .1]';
% 
% psi{1} = [.1 .7 .2; 
%           .3 .6 .1];
% psi{2} = [.2 .8; 
%           .6 .4];
% psi{3} = [.1 .2 .65 .05;
%           .1 .7 .05 .15];
%
% L = loglik(seq, E, phi, psi);
%
% Author: geopapa
% $ Date: 2013/10/15 22:25:23 $

    N  = length(seq);
    Ne = size(E,1);
    L  = 0;
    if iscell(phi)
        for i = 1:N  % Sum-up the node potentials
            L = L + log(phi{i}(seq(i)));
        end
    else
        for i = 1:N  % Sum-up the node potentials
            L = L + log(phi(seq(i)));
        end
    end
    
    if iscell(psi)
        for e = 1:Ne  % Sum-up the pairwise potentials
            L = L + log(psi{e}(seq(E(e,1)),seq(E(e,2))));
        end
    else
        for e = 1:Ne  % Sum-up the pairwise potentials
            L = L + log(psi(seq(E(e,1)),seq(E(e,2))));
        end
    end
end