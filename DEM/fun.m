function err = fun(x,xi,y)
covar.model = 'k-bessel';
covar.alpha = x(1);
covar.var = x(2);
covar.range = [x(3) x(3)];
covar = covarIni(covar);

err = nansum( ( 1-covar.g(xi*covar.cx(1)) - y' ).^2 );

end