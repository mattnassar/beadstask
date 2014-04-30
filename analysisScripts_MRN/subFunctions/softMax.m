function [pChoice]=softMax(values, invT)

values=values-(max(values));

pChoice=exp(values.*invT)./nansum(exp(values.*invT));
