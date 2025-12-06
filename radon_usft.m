function R=radon_usft(f,theta)
N=size(f,1);
Ntheta=size(theta,2);

%circ support
[xx,yy]=ndgrid((-N/2:N/2-1),(-N/2:N/2-1));
circ=(2/N*xx).^2+(2/N*yy).^2<1-1/N;
f=f.*circ;

%polar grid
rho=(-N/2:N/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);

%border control
x(find(x>=0.5))=0.5-1e-8;
y(find(y>=0.5))=0.5-1e-8;

% fourier transform on unequally-spaced grid
[xeq,yeq] = meshgrid(-N/2:N/2-1);
F=zeros(size(x));
for j0=1:Ntheta
    for j1=1:N
        F(j1,j0) = F(j1,j0)+sum(sum(f .* exp(-2*pi*1j*(yeq*y(j1,j0)+xeq*x(j1,j0)))));        
    end
end

%IFFT1D
R=fftshift(ifft(ifftshift(F)));
