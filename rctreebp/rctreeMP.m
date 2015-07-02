% rctreeMP(A,F,K) : build an RC-tree structure for max-product from a factor graph
%   A: sparse adjacency matrix, Nv x Nf  (#vars by #factors)
%   F: cell array of factors, 1 x Nf
%      each factor F{i} is a matrix D1xD2xD3x..., the dimensions of the adjacent variables in lexicographic order
%      so that if F{1} is adj to v1, v4, and v2, F{1} is |v1| x |v2| x |v4|
%   K: scalar integer controlling the amount of memory allocated for RCtree growth
%      K too small can force memory re-allocation; K too large wastes memory.
%
function T=rctreeMP(A,F,K)

if (nargin < 3) K=4; end;
N=size(A,1); M=size(A,2);  NN=N+M;  Nrounds=ceil(K*log2(NN));

%%%%%%% DATA STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T.First = zeros(1,Nrounds,'uint32'); % Nrounds = N?
T.Latest = zeros(1,NN,'uint32');
T.MemPtr = uint32(1);

T.Next=zeros(1,K*NN,'uint32');  T.Prev=zeros(1,K*NN,'uint32');
T.Newer=zeros(1,K*NN,'uint32'); T.Older=zeros(1,K*NN,'uint32');
T.Id=zeros(1,K*NN,'uint32');    T.Roll=zeros(1,K*NN,'uint8');
T.Adj=cell(1,K*NN);
T.in=zeros(1,NN,'uint32'); T.out=zeros(1,NN,'uint32');
T.Var=cell(1,NN); T.F=cell(1,NN); T.Dim=zeros(1,N,'uint32'); 
T.CF=cell(1,NN); T.CFVar=cell(1,NN);
T.Scars = cell(1,K*NN);
T.K=K;				% current multiplicative factor
T.roots=uint32([]);		% roots of RCtree forest
T.map = zeros(1,N,'uint32');	% current map configuration of vars
T.diff = zeros(1,N,'uint32');
T.changed=zeros(1,N,'uint32');
T.changedCount=uint32(0);
DEBUG = 0;			% no informational / debugging messages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE using other adj matrix form...
%  This can take a while, but we don't count it against the algorithm since it's
%  Matlab's slowness of converting a sparse matrix into a pointer-based structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=N+1:N+M,
  T.Next(i)=i+1; T.Prev(i)=i-1; T.Newer(i)=0; T.Older(i)=0; T.Id(i)=i; T.Roll(i)=rand>.5;
  T.Adj{i}=uint32(full(find(A(:,i-N)))'); T.Var{i}=T.Adj{i}; T.F{i}=F{i-N};
%  for j=T.Adj{i}, T.Adj{j}=union_sorted(T.Adj{j},i); end;
  tmp=size(F{i-N}); T.Dim(T.Var{i})=tmp(1:length(T.Var{i}));
  if (length(T.Adj{i})<2) T.Roll(i)=T.Roll(i)+2; end;  % mark leaves/roots as +2
end;
for i=1:N,
  T.Next(i)=i+1; T.Prev(i)=i-1; T.Newer(i)=0; T.Older(i)=0; T.Id(i)=i; T.Roll(i)=rand>.5;
  T.Adj{i} = uint32(full(find(A(i,:)))+N); 
  T.Var{i}=uint32(i); T.F{i}=ones(T.Dim(i),1);
  if (length(T.Adj{i})<2) T.Roll(i)=T.Roll(i)+2; end;  % mark leaves/roots as +2
end;
T.First(1)=1; T.Latest=uint32(1:N+M); 
T.MemPtr=NN+1; T.Next(NN)=0; T.Next(NN+1:K*NN-1)=NN+2:K*NN;
% % % tic;   % I HAVE COMMENTED THIS OUT CAUSE IT IS NOT FAIR TO NOT COUNT
             % IT IN CALCULATIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin constructing the RCTree by rake & compress operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
round=1; while(T.First(round)~=0),

  if (DEBUG)	% For debugging: print out all the current unraked nodes
    ptr=T.First(round); while (ptr~=0), fprintf('%d ',T.Id(ptr)); ptr=T.Next(ptr); end; fprintf('\n');
  end;

  %%%%%%%%%%%%%%%%%%% DECIDE WHICH NODES TO KEEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ptr=T.First(round); while (ptr~=0), % if we're going to keep it,
    if ( length(T.Adj{ptr})>2 || ...  % keep: all 3+deg; some 2-deg; one of leaf-leaf pairs 
        (length(T.Adj{ptr})==2 && (T.Roll(ptr)==0 || any(T.Roll(T.Adj{ptr}))) ) || ...
	(length(T.Adj{ptr})==1 && length(T.Adj{T.Adj{ptr}})==1 && T.Adj{ptr}<ptr )  ),      
      % copy it into the next set

      %%%%%%%%%%%%%%%%%%% ALLOCATE MEMORY FOR KEPT NODES %%%%%%%%%%%%%%%%%%%%%
      if (T.MemPtr==0)                        % OUT OF MEMORY?
        if (DEBUG) fprintf('Out of memory...\n'); end;
        K=T.K; T.K=K+1; NN=length(T.in); KN=K*NN; KKN=KN+NN;
        T.Next(KKN)=0; T.Prev(KKN)=0;         % If so, extend everything by NN
        T.Newer(KKN)=0; T.Older(KKN)=0;
        T.Id(KKN)=0; T.Roll(KKN)=0;
        T.Adj{KKN}=[];
        T.Next(KN:KKN-1)=KN+1:KKN; T.MemPtr=KN;
      end;
      new=T.MemPtr; T.MemPtr=T.Next(T.MemPtr);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      T.Next(new)=T.First(round+1); T.First(round+1)=new;
      if (T.Next(new)~=0) T.Prev(T.Next(new))=new; end;
      T.Latest(T.Id(ptr))=new; T.Newer(ptr)=new; T.Older(new)=ptr; T.Newer(new)=0;
      T.Roll(new)=rand>.5; T.Id(new)=T.Id(ptr);
      T.Scars{new}=T.Scars{ptr};
    else T.Newer(ptr)=0;
    end;
    ptr=T.Next(ptr);
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%% FOR EACH NEWLY DELETED NODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ptr=T.First(round); while (ptr~=0), if (T.Newer(ptr)==0),
      % deleting the node -- get our scars, remove from circulation, add self
      nbrs=T.Adj{ptr}; v=T.Id(ptr); cfs=T.Scars{ptr}; % 
      %for i=T.Id(nbrs), T.Scars{i}=union(setdiff(T.Scars{i},cfs),ptr); end;
      for i=nbrs, T.Scars{T.Newer(i)}=sUnion32(sDiff32(sort(T.Scars{T.Newer(i)}),sort(cfs)),ptr); end;
      %%%%%%%%%%%%%%%%%%% COMPUTE NEW CLUSTER FUNCTION CF %%%%%%%%%%%%%%%%%%%%%
      vars=T.Var{v};    % get variables and their dimensions...
      for i=T.Id(cfs), vars = sUnion32(vars, T.CFVar{i}); end;
      dims=T.Dim(vars);
      CF = multCF(1,vars,dims,T.F{v},T.Var{v});
      for i=T.Id(cfs), CF = multCF(CF,vars,dims,T.CF{i},T.CFVar{i}); end;

      if (DEBUG)	% For debugging: print out the current operation
        switch (length(nbrs))
        case 0, fprintf('Rooting node %d\n',v);
        case 1, fprintf('Raking node %d [%d]\n',v,T.Id(nbrs(1)));
        case 2, fprintf('Compressing node %d [%d %d]\n',v,T.Id(nbrs(1)),T.Id(nbrs(2)));
        otherwise, fprintf('Improper removal of node %d\n',v);
        end;	   
      end;
      
      switch (length(nbrs))
      case 0, CF=max(CF(:)); s1=[]; T.roots=sUnion32(T.roots,v); keep=[];
      case 1, T.in(v)=T.Id(nbrs); keep=T.Var{T.Id(nbrs)};
      case 2, T.in(v)=T.Id(nbrs(1)); T.out(v)=T.Id(nbrs(2));
	      keep=sUnion32(T.Var{T.Id(nbrs(1))},T.Var{T.Id(nbrs(2))});
      end;
      hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
      if (~isscalar(vars)), CF=permute(CF,[i1,i2]); end;
      for i=length(vars):-1:length(s1)+1, CF=max(CF,[],i); end;
      T.CF{v} = CF; T.CFVar{v}=s1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end; ptr=T.Next(ptr); end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%% FOR EACH REMAINING NODE FIX ADJACENCY ETC %%%%%%%%%%%%%%%
  round = round+1;
  ptr=T.First(round); while (ptr~=0)
    % Fix the adjacency structure of the new round
    old=T.Older(ptr); AdjO=T.Adj{old}; AdjN=T.Newer(AdjO);
    for i=find(AdjN==0),
      tmp = T.Newer(T.Adj{AdjO(i)}); 
      if (tmp(1)~=ptr), AdjN(i)=tmp(1); elseif (length(tmp)>1) AdjN(i)=tmp(2); end;
    end;
    AdjN=AdjN(find(AdjN~=0)); T.Adj{ptr}=AdjN;
    if (length(AdjN)<2) T.Roll(ptr)=T.Roll(ptr)+2; end;  % mark as leaf or root...
    ptr=T.Next(ptr);
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;

% FIX UP IN/OUT DIRECTIONS AND COMPUTE MAP
D=zeros(1,NN,'uint32');				% vector of "visited" tags for postprocessing
for v=1:NN, postprocess(T,v,D); end;		% compute in&out relationships & initial MAP config
T.diff = zeros(1,N,'uint32');			% reset MAP change detection variables
T.changed=zeros(1,N,'uint32');
T.changedCount=uint32(0);




function postprocess(T,v,D)

  if (v > 0 && D(v) == 0)			% no map estimate; not yet done
    u=T.in(v); w=T.out(v);
    postprocess(T,u,D); postprocess(T,w,D);	% recurse on parents
    if (w ~= 0) 
     if ( (T.in(u)==w) || (T.out(w)==u) ),  
      modU32Array(T.in,v,w); modU32Array(T.out,v,u); 	%T.in(v)=w; T.out(v)=u;
     else 
      modU32Array(T.in,v,u); modU32Array(T.out,v,w); 	%T.in(v)=u; T.out(v)=w;
     end;
    end;
    updateMAP(T,v);				% update the MAP configuration
    %[vals,vars] = maxmarginalC(T,v);
    %vals=uint32(vals); vars=uint32(vars);
    %for i=1:length(vars), modU32Array(T.map,uint32(vars(i)),uint32(vals(i))); end;
    modU32Array(D,v,uint32(1));			% mark this node as complete
  end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute "orientations" of nodes in the RC Tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function orient(v)
  global RCTREE_IN RCTREE_OUT RCTREE_TODO
  if (RCTREE_OUT(v)~=0)                       % if v was compressed, in & out store nbrs
    u=RCTREE_IN(v); w=RCTREE_OUT(v);          %  but not the correct orientation
    if (RCTREE_TODO(u)), orient(u); end;      % if nbrs not yet oriented, do so
    if (RCTREE_TODO(w)), orient(w); end;      %  then orient v to its nbrs
    if ((RCTREE_IN(u)==w)||(RCTREE_OUT(w)==u)), RCTREE_IN(v)=w; RCTREE_OUT(v)=u; 
    else RCTREE_IN(v)=u; RCTREE_OUT(v)=w;     %  using their orient. (v is between them)
    end;
  end; 
  RCTREE_TODO(v)=0;			      % mark as done


