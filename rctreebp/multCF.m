% CF=multCF(CF,all,dims,cf,vs) -- internal function to multiply cluster f'ns
%  given a cf over "all" variables with dimensions dim, multiply it by a new
%  factor defined over only the (subset of) variables vs.
%
function CF=multCF(CF,allvar,dims,cfi,cfVi)
  if (isempty(cfVi)), return; end;		% cfi not a function?
  if (length(allvar)==1), CF = CF .* cfi;	% or everything is univariate? 
  else
   hits = sMember32(allvar,cfVi);
%  [incl,iidx]=intersect(allvar,cfVi); 		% which variables is cfi over
%  [excl,eidx]=setdiff(allvar,incl);		%   and which is it not?
%  drep=dims; dresh=dims; dresh(eidx)=1; drep(iidx)=1;  % reshape it to the
  drep=double(dims); dresh=drep; dresh(~hits)=1; drep(hits)=1;  % reshape it to the
  CF = CF .* repmat(reshape(cfi,dresh), drep);  %   right dimension & copy
  end;						%   values across unused dims

