function f = evalgamma(params, X)

f = params(1).*gampdf(X,params(2),params(3));