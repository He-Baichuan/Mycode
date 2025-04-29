function S = contemp(X,name,m)
    idx =  ModelUtils.id(name,m,0);
    S = X(idx:m.nx:end);
end