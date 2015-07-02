function M = rmq(A)
% function M = rmq(A)
% 
% This function gives the RMQ structure of an array A <http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=lowestCommonAncestor LCA to RMQ>
% It builds the RMQ structure of an array A, which allows for the recovery 
% of the index of the minimum element of any subarray of A in constant time.
%
% INPUT
%  A:  1xN   array     -  Input array      
% 
% 
% OUTPUT
%  M:  NxlogN  matrix  -  RMQ structure
%                         M(i,j): Absolute index of minimum element in 
%                         A[i...i+2^k-1], where k = floor(log2(j-i+1)). 
%                         In other words, it gives the absolute index in
%                         the interval of length 2^k, which starts at i
% 
% 
% EXAMPLES
%
% A = [-1  5  0  3  -6  4  2  1  0  -1];
% N = length(A);
% M = rmq(A);
% % Determine the index of the minimum element in A[i...j] in O(1) time
% for i = 1:N
%     for j = i:N
%         k = floor(log2(j - i + 1));
%         if A(M(i,k+1)) <= A(M(j-2^k+1,k+1))
%             idx = M(i,k+1);
%         else
%             idx = M(j-2^k+1,k+1);
%         end
%     end
% end
% 
% The above example shows how we can build the RMQ structure M and recover
% minimum index queries in constant O(1) time.
% 
% Complexity of building RMQ structure: time: O(N*logN), space: O(N*logN)
% Complexity of finding the index of min element in A[i...j]: time: O(1)
% 
% Author: geopapa
% $ Date: 2014/01/03 12:01:32 $    
    
    N = length(A);
    L = ceil(log2(N))+1;
    M = zeros(N,L);
 
    % Initialize M for the intervals with length 1
    M(:,1) = (1:N)';
    
    % Compute values of M dynamically
    for j = 2:L
        k = bitshift(1,(j-1)-1);  % j-1 is the previous column and -1 is due to the fact that (j-1)^th col corresponds to 2^(j-1-1) range
        for i = 1:N
            if A(M(i,j-1)) <= A(M(min(i+k,N),j-1))
                M(i,j) = M(i,j-1);
            else
                M(i,j) = M(min(i+k,N),j-1);
            end
        end
    end
end

function idx = query(i, j, A, M)
    k = floor(log2(j - i + 1));
    if A(M(i,k+1)) <= A(M(j-2^k+1,k+1))
        idx = M(i,k+1);
    else
        idx = M(j-2^k+1,k+1);
    end
end