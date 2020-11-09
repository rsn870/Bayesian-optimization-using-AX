function [xxnew] = proposal(xx,ntot,ep)
n_xx = length(xx);
n_replace = ep*n_xx;
indx_replace = randperm(n_xx,n_replace);
val_replace = randperm(ntot,n_replace);
xxnew = xx;
xxnew(indx_replace)=val_replace;
xxnew = unique(xxnew);

n_xxnew = length(xxnew);


while n_xxnew<n_xx
    ndiff = n_xx-n_xxnew;
    val_add = randperm(ntot,ndiff);
    xxnew = [xxnew,val_add];
    xxnew = unique(xxnew);
    n_xxnew = length(xxnew);
end
    