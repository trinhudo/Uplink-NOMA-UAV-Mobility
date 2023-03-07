function output...
    = I1(K0,chi0,Omg0,mu0,p0,q0,...
         K1,chi1,Omg1,mu1,p1,q1)
    
    result = 0;
    for k0 = 1:(K0+1)
    for k1 = 1:(K1+1)
       kap0= mu0(k0);
    for i0 = 0:(kap0-1)
    for j0 = 0:i0
       l0 = i0-j0;
       %
       Lambda0 = Omg0(k0);
       Lambda1 = (1/Omg1(k1)+p0/Lambda0)^(-1);
       kap1= mu1(k1) + l0;
       Xi0 = @(k) chi0(k)/factorial(mu0(k)-1) * Omg0(k)^(-mu0(k));
       Xi1 = @(k) chi1(k)/factorial(mu1(k)-1) * Omg1(k)^(-mu1(k));
       %
       result = result + Xi0(k0) * Xi1(k1) * Lambda0^(-i0)...
           * q0^j0/factorial(j0) * p0^l0/factorial(l0) * Lambda1^kap1...
           * gammainc((p1+q1)/Lambda1,kap1,'upper')*gamma(kap1)...
           * factorial(mu0(k0)-1) * exp( -q0/Lambda0 )*Lambda0^kap0;
       %
    end
    end
    end
    end
    %
    output = result;
end