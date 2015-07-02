% [p,v] = belief(T,i) -- compute beliefs at some node i in the RCTree
%   p is the normalized belief (marginal) over variables v.  If i is a
%   variable node, p is a column vector and v=i; if i is a factor node,
%   p is a matrix (joint distribution) and v a list of involved vars.
%
function [p,var]=marginal(T,v)
  v = uint32(v);				% convert to uint just in case
  orient(T,v);					% orient the nodes above v
  %fprintf('Computing belief %d\n',v);

  % recursively find the messages from our parent nodes, following "out" first.
  [msgout,vout,lastOutPar] = message(T,T.out(v),v, 0,v);
  [msgin,vin] = message(T,T.in(v),v, 0,lastOutPar);

  % get list of all involved variables in our marginal computation
  vars=sUnion32(vout,vin); vars=sUnion32(vars,T.Var{v}); cfs=T.Id(T.Scars{T.Latest(v)});
  for i=cfs, vars=sUnion32(vars, T.CFVar{i}); end; dims=T.Dim(vars);

  % compute product of all cluster f'ns, incoming msgs, etc:
  CF = multCF(1,vars,dims,msgout,vout);  
  CF = multCF(CF,vars,dims,msgin,vin);
  CF = multCF(CF,vars,dims,T.F{v},T.Var{v});
  for i=1:length(cfs)
    CF = multCF(CF,vars,dims,T.CF{cfs(i)},T.CFVar{cfs(i)});
  end;
  keep = T.Var{v};  				% marginalize over "other" variables
  hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
  %[s1,i1]=intersect(vars,keep); [s2,i2]=setdiff(vars,keep);
  if (~isscalar(vars)) CF=permute(CF,[i1,i2]); end;
  for i=length(vars):-1:length(s1)+1, CF=sum(CF,i); end;
  p=CF/sum(CF(:)); var=s1;     % set function outputs



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Compute/update the orientation (in/out) of the nodes above v
%
function orient(T,v)
  if (T.in(v)==0) return; end;			% if we're a root, we're done
  u=T.in(v); w=T.out(v);			% else, orient our parent (lowest of u,w)
  if (w==0 || T.in(u)==w || T.out(u)==w) orient(T,u); else orient(T,w); end;
  if (w~=0 && (T.in(u)==w || T.out(w)==u)), 	% then find our orientation given our parents
    modU32Array(T.in,v,w); modU32Array(T.out,v,u);  %T.in(v)=w; T.out(v)=u; 
  else 
    modU32Array(T.in,v,u); modU32Array(T.out,v,w);  %T.in(v)=u; T.out(v)=w; 
  end;
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Compute (downward) RC-tree message m_{u->v}
%   "followOuts" = binary indicator: do we need the message from u's out (if it exists)?
%   "lastOutPar" = the node whom our last "out" message was raked into.  From this point on,
%                   we'll need to look for "out" messages again.
%
function [msg,var,lastOutPar]=message(T,u,v,followOuts,lastOutPar)
  if (u==0) 					% terminate recursion: 
    msg=1; var=[]; lastOutPar=T.in(v); 		%  if we were following "out" paths, v was
    return; 					%  raked into T.in(v), so return that as well.
  end;
  %fprintf('Computing message %d -> %d\n',u,v);
  % follow recursion to get incoming message(s) -- note: must always check "out" paths before "in" paths
  followOuts = (followOuts || u==lastOutPar); 	% do we need to follow an out path, if it exists?
  if (T.out(v)==u) 				% if we're following an out path right now,
    [msg,var,lastOutPar]=message(T,T.out(u),u ,0,u);	% keep following it.
  else						% if we're on an in-path and looking, check out also
    if (followOuts) [msgout,vout,lastOutPar]=message(T,T.out(u),u, followOuts,lastOutPar);
                    [msgin,vin]=message(T,T.in(u),u, 0,lastOutPar);
                    var=sUnion32(vout,vin);	%   need the product of in & out messages
                    msg=multCF(1,var,T.Dim(var),msgout,vout); 
                    msg=multCF(msg,var,T.Dim(var),msgin,vin);
    else 					%   if we're not looking for out paths, just stick
      [msg,var]=message(T,T.in(u),u, followOuts,lastOutPar);	% to the in-path.
    end;
  end;
  tmp=[0,v]; 					% find the direction of v from u (u's child)
  while (tmp(2)~=u), 				%   we'll leave that cluster f'n out of the msg computation
    tmp(1)=tmp(2);  uu=T.in(tmp(1)); ww=T.out(tmp(1));
    if ((ww==0)||(T.in(uu)==ww)||(T.out(uu)==ww)) tmp(2)=uu; else tmp(2)=ww; end;
  end; chEx=tmp(1);				%   (leave out chEx from product for msg)

  % get list of all involved variables:
  vars=sUnion32(var,T.Var{u}); 			% get "cfs", the clusters below u (not including chEx)
  cfs=sDiff32( sort(T.Id(T.Scars{T.Latest(u)})) , chEx );
  for i=cfs, vars=sUnion32(vars, T.CFVar{i}); end;	% and find all the used variables.
  dims=T.Dim(vars);

  % compute product of all cluster f'ns, etc.:
  CF=multCF(1,vars,dims,msg,var); CF=multCF(CF,vars,dims,T.F{u},T.Var{u});
  for i=1:length(cfs)
    CF = multCF(CF,vars,dims,T.CF{cfs(i)},T.CFVar{cfs(i)});
  end;
  keep = T.CFVar{v};   				% marginalize over "other" variables
  hits=sMember32(vars,keep); i1=find(hits); s1=vars(i1); i2=find(~hits);
  %[s1,i1]=intersect(vars,keep); [s2,i2]=setdiff(vars,keep);
  if (~isscalar(vars)) CF=permute(CF,[i1,i2]); end;
  for i=length(vars):-1:length(s1)+1, CF=sum(CF,i); end;
  msg=CF/sum(CF(:)); var=s1;    		% set function outputs

  %fprintf('Done message %d -> %d\n',u,v); msg,


