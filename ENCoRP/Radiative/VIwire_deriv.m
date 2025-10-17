function out = VIwire_deriv(x,y,z,wpw,f,whichxyz,VorC)
    out = zeros(length(x),1);
    h = 1e-3; % h is currently hardcoded into VIwire
    for i = 1:length(x)
        l = VIwire(x(i),y(i),z(i),wpw(i),-1,f(i),VorC(i));
        r = VIwire(x(i),y(i),z(i),wpw(i),1,f(i),VorC(i));
        dR = (r-l)/(2*h);
        out(i) = dR(whichxyz(i));
    end
    out(isnan(out)) = 0;
    out(isinf(out)) = 0;
end