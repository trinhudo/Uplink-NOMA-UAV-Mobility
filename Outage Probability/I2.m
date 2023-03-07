function output...
    = I2(K0,chi0,Omg0,mu0,p0,q0,...
         K1,chi1,Omg1,mu1,p1,q1,...
         K2,chi2,Omg2,mu2,p2,q2)
    
%     result = 0;
%     for k0 = 1:(K0+1)
%     for k1 = 1:(K1+1)
%     for k2 = 1:(K2+1)
%        kap0= mu0(k0);
%     for i0 = 0:(kap0-1)
%     for j0 = 0:i0
%        l0 = i0-j0;
%        kap1= mu1(k1) + l0;
%     for i1 = 0:(kap1-1)
%     for j1 = 0:i1
%        l1 = i1-j1;
%        kap2= mu2(k2) + l1;
       %
       Lambda0 = Omg0(k0);
       Lambda1 = (1/Omg1(k1)+p0/Lambda0)^(-1);
       Lambda2 = (1/Omg2(k2)+p1/Lambda1)^(-1);
       %
       Xi0 = @(k) chi0(k)/factorial(mu0(k)-1) * Omg0(k)^(-mu0(k));
       Xi1 = @(k) chi1(k)/factorial(mu1(k)-1) * Omg1(k)^(-mu1(k));
       Xi2 = @(k) chi2(k)/factorial(mu2(k)-1) * Omg2(k)^(-mu2(k));
       %
       result = Xi0(k0) * Xi1(k1) * Xi2(k2) * Lambda0^(-i0)...
           * q0^j0/factorial(j0) * p0^l0/factorial(l0) * Lambda1^(-i1)...
           * q1^j1/factorial(j1) * p1^l1/factorial(l1) * Lambda2^(kap2)...
           * gammainc((p2+q2)/Lambda2,kap2,'upper')*gamma(kap2)...
           * factorial(kap0-1) * exp( -q0/Lambda0 )*Lambda0^kap0...
           * factorial(kap1-1) * exp( -q1/Lambda1 )*Lambda1^kap1;
%        result = result + Xi0(k0) * Xi1(k1) * Xi2(k2) * Lambda0^(-i0)...
%            * q0^j0/factorial(j0) * p0^l0/factorial(l0) * Lambda1^(-i1)...
%            * q1^j1/factorial(j1) * p1^l1/factorial(l1) * Lambda2^(kap2)...
%            * gammainc((p2+q2)/Lambda2,kap2,'upper')*gamma(kap2)...
%            * factorial(kap0-1) * exp( -q0/Lambda0 )*Lambda0^kap0...
%            * factorial(kap1-1) * exp( -q1/Lambda1 )*Lambda1^kap1;
       %
%     end
%     end
%     end
%     end
%     end
%     end
%     end
    %
    output = result;
end