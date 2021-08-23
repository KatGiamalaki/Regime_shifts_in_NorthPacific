function [theta]=extremal_Sueveges_events(Y,p)
u=quantile(Y, p);
q=1-p;
Li=find(Y>u);
Ti=diff(Li);
Si=Ti-1;
Nc=length(find(Si>0));
N=length(Ti);

theta=(sum(q.*Si)+N+Nc- sqrt( (sum(q.*Si) +N+Nc).^2  -8*Nc*sum(q.*Si))  )./(2*sum(q.*Si));
end