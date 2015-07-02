function subs=ind2sub2(dims,indx)
  subs=ones(size(dims)); indx=double(indx-1); dims=double(dims);
  subs(1)=1+mod(indx,dims(1));
  tmp=dims(1);
  for i=2:length(dims)
    subs(i)=floor(indx/tmp)+1;
    tmp=tmp*dims(i);
  end;




