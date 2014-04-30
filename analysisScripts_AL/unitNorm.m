function [y]=unitNorm(x)

y=(x-min(x))./max(x-min(x));