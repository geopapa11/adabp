function T=rctreeSP(A,F,K)

if (nargin < 3) K=4; end;
N=size(A,1); M=size(A,2);  NN=N+M;  Nrounds=ceil(K*log2(NN));
DEBUG=0;
%% DATA STRUCTURE
T.First = zeros(1,Nrounds,'uint32'); 		% beginning of linked list for each round of tree contraction
T.Latest = zeros(1,NN,'uint32');		% last copy of node i before it is removed during contraction
T.MemPtr = uint32(1);				% next free memory location in heap

T.Next=zeros(1,K*NN,'uint32');  T.Prev=zeros(1,K*NN,'uint32');	% linked list data structure for each round
T.Newer=zeros(1,K*NN,'uint32'); T.Older=zeros(1,K*NN,'uint32');	% node i's connections to prior & subseq rounds
T.Id=zeros(1,K*NN,'uint32');    T.Roll=zeros(1,K*NN,'uint8');	% the ID # and random coin flip for each node
T.Adj=cell(1,K*NN);						% (partially contracted) adjacency structure of each node 
T.in=zeros(1,NN,'uint32'); T.out=zeros(1,NN,'uint32');		% in & out neighbors (neighbors on removal)
T.Var=cell(1,NN); T.F=cell(1,NN); T.Dim=zeros(1,N,'uint32'); 	% functions associated with each *node* of G
T.CF=cell(1,NN); T.CFVar=cell(1,NN);				% functions associated with each *cluster* at removal
T.Scars = cell(1,K*NN);						% removed neighbors during contraction
T.K=K;								% memory overhead factor
T.roots=uint32([]);						% list of roots for the RCtree forest


%tic;
% INITIALIZE using other adj matrix form... (this can take a while, due to Matlab's sparse matrix implementation)
for i=N+1:N+M,				% init factors
  T.Next(i)=i+1; T.Prev(i)=i-1; T.Newer(i)=0; T.Older(i)=0; T.Id(i)=i; T.Roll(i)=rand>.5;
  T.Adj{i}=uint32(full(find(A(:,i-N)))'); T.Var{i}=T.Adj{i}; T.F{i}=F{i-N};
%  for j=T.Adj{i}, T.Adj{j}=union_sorted(T.Adj{j},i); end;
  tmp=size(F{i-N}); T.Dim(T.Var{i})=tmp(1:length(T.Var{i}));
  if (length(T.Adj{i})<2) T.Roll(i)=T.Roll(i)+2; end;  % mark leaves/roots as +2
end;
for i=1:N,				% init variables
  T.Next(i)=i+1; T.Prev(i)=i-1; T.Newer(i)=0; T.Older(i)=0; T.Id(i)=i; T.Roll(i)=rand>.5;
  T.Adj{i} = uint32(full(find(A(i,:)))+N); 
  T.Var{i}=uint32(i); T.F{i}=ones(T.Dim(i),1);
  if (length(T.Adj{i})<2) T.Roll(i)=T.Roll(i)+2; end;  % mark leaves/roots as +2
end;
T.First(1)=1; T.Latest=uint32(1:N+M); 	% init first round of linked list and free memory
T.MemPtr=NN+1; T.Next(NN)=0; T.Next(NN+1:K*NN-1)=NN+2:K*NN;
% % % tic;   % I HAVE COMMENTED THIS OUT CAUSE IT IS NOT FAIR TO NOT COUNT
             % IT IN CALCULATIONS

round=1; while(T.First(round)~=0),	% Each round of tree contraction:
  if (DEBUG) ptr=T.First(round); while (ptr~=0), fprintf('%d ',T.Id(ptr)); ptr=T.Next(ptr); end; fprintf('\n'); end;
  ptr=T.First(round); while (ptr~=0), 	% if we're going to keep it,
    if ( length(T.Adj{ptr})>2 || ... 	% (keep: all 3+deg, some 2-deg, and one of leaf-leaf pairs)
        (length(T.Adj{ptr})==2 && (T.Roll(ptr)==0 || any(T.Roll(T.Adj{ptr}))) ) || ...
	(length(T.Adj{ptr})==1 && length(T.Adj{T.Adj{ptr}})==1 && T.Adj{ptr}<ptr )  ),      
      					% copy it into the next set
      if (T.MemPtr==0)                        % OUT OF MEMORY?
        %fprintf('Out of memory...\n'); 
        K=T.K; T.K=K+1; NN=length(T.in); KN=K*NN; KKN=KN+NN;	
        T.Next(KKN)=0; T.Prev(KKN)=0;		% If so, extend everything by NN 
        T.Newer(KKN)=0; T.Older(KKN)=0; 
        T.Id(KKN)=0; T.Roll(KKN)=0; 
        T.Adj{KKN}=[];
        T.Next(KN:KKN-1)=KN+1:KKN; T.MemPtr=KN; 
        end;            
      new=T.MemPtr; T.MemPtr=T.Next(T.MemPtr);
      %fprintf('Allocating memory: node %d\n',new);

      T.Next(new)=T.First(round+1); T.First(round+1)=new;
      if (T.Next(new)~=0) T.Prev(T.Next(new))=new; end;
      T.Latest(T.Id(ptr))=new; T.Newer(ptr)=new; T.Older(new)=ptr; T.Newer(new)=0;
      T.Roll(new)=rand>.5; T.Id(new)=T.Id(ptr);
      T.Scars{new}=T.Scars{ptr};
    else T.Newer(ptr)=0;		% if we removed it, mark it for removal
    end;
    ptr=T.Next(ptr);
  end;

  ptr=T.First(round); while (ptr~=0), if (T.Newer(ptr)==0),
      % deleting the node -- get our scars, remove from circulation, add self
      nbrs=T.Adj{ptr}; v=T.Id(ptr); cfs=T.Scars{ptr}; % 
      for i=nbrs, T.Scars{T.Newer(i)}=sUnion32(sDiff32(sort(T.Scars{T.Newer(i)}),sort(cfs)),ptr); end; % !!! sort needed?
      %%%%%%%%%%%%%%%%%%% COMPUTE NEW CF FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vars=T.Var{v};    % get variables and their dimensions...
      for i=T.Id(cfs), vars = sUnion32(vars, T.CFVar{i}); end;
      dims=T.Dim(vars);
      CF = multCF(1,vars,dims,T.F{v},T.Var{v});
      for i=T.Id(cfs), CF = multCF(CF,vars,dims,T.CF{i},T.CFVar{i}); end;

      if (DEBUG), switch (length(nbrs))
        case 0, fprintf('Rooting node %d\n',v);
        case 1, fprintf('Raking node %d [%d]\n',v,T.Id(nbrs(1)));
        case 2, fprintf('Compressing node %d [%d %d]\n',v,T.Id(nbrs(1)),T.Id(nbrs(2)));
        otherwise, fprintf('Improper removal of node %d\n',v);
      end; end;
	      
      switch (length(nbrs))
      case 0, CF=sum(CF(:)); s1=[]; T.roots=sUnion32(T.roots,v); keep=[];	% created a new root
      case 1, T.in(v)=T.Id(nbrs); keep=T.Var{T.Id(nbrs)};			% raked v 
      case 2, T.in(v)=T.Id(nbrs(1)); T.out(v)=T.Id(nbrs(2));			% compressed v
	      keep=sUnion32(T.Var{T.Id(nbrs(1))},T.Var{T.Id(nbrs(2))});
      end;
      hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);	% marginalize out the
      if (~isscalar(vars)), CF=permute(CF,[i1,i2]); end;			% cluster function on 
      for i=length(vars):-1:length(s1)+1, CF=sum(CF,i); end;			% "other" variables
      T.CF{v} = CF; T.CFVar{v}=s1;						% save function CF(s1)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end; ptr=T.Next(ptr); end;

  round = round+1;						% get ready for the next round:
  ptr=T.First(round); while (ptr~=0)				% fix up the adjacency structure
    old=T.Older(ptr); AdjO=T.Adj{old}; AdjN=T.Newer(AdjO);	% get pointers to the old & new versions of our nbrs
    for i=find(AdjN==0),					% any of them that were removed:
      tmp = T.Newer(T.Adj{AdjO(i)}); 				%  find *their* old neighbors in this round
      if (tmp(1)~=ptr), AdjN(i)=tmp(1); 			%  removed nodes had at most two neighbors, so find
      elseif (length(tmp)>1) AdjN(i)=tmp(2); end;		%  the one that's not "ptr" and connect to it.
    end;
    AdjN=AdjN(find(AdjN~=0)); T.Adj{ptr}=AdjN;			% save the new adjacency list ("auto" sorted)
    if (length(AdjN)<2) T.Roll(ptr)=T.Roll(ptr)+2; end;  	% mark the node as being a leaf or root...
    ptr=T.Next(ptr);
  end;
end;

%% FIX UP IN/OUT DIRECTIONS...					% this part is no longer done here
%global RCTREE_IN RCTREE_OUT RCTREE_TODO			% we now compute direction in the query function
%RCTREE_IN=T.in; RCTREE_OUT=T.out; RCTREE_TODO=T.out~=0;
%for v=1:NN, if (RCTREE_TODO(v)), orient(v); end; end;
%T.in=RCTREE_IN; T.out=RCTREE_OUT;
%clear global RCTREE_IN RCTREE_OUT RCTREE_TODO




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



