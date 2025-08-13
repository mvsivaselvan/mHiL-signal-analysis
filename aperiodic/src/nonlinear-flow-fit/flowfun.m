function q = flowfun(xv)

k1 = 1.18e4;
k2 = 3.2e4;
delt = 0.25;

q = zeros(size(xv));

q(xv<=delt) = k1*xv(xv<=delt);
q(xv>delt) = k1*delt + k2*(xv(xv>delt) - delt);
q(xv<-delt) = -k1*delt + k2*(xv(xv<-delt) + delt);
