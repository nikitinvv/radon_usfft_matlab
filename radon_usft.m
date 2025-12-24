function R=radon_usft(f,theta,os)
N=size(f,1);
Ntheta=size(theta,2);

%polar grid
rho=(-N*os/2:N*os/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);

% fourier transform on unequally-spaced grid
[xeq,yeq] = meshgrid(-N/2:N/2-1);
F=zeros(size(x));
for j0=1:Ntheta
    for j1=1:N*os
        F(j1,j0) = F(j1,j0)+sum(sum(f .* exp(-2*pi*1j*(yeq*y(j1,j0)+xeq*x(j1,j0)))));        
    end
end

%IFFT1D
R=fftshift(ifft(ifftshift(F)));

% bin
if (os==2)
    R = R(1:os:end,:)+R(2:os:end,:);
end