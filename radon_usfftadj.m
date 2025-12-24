function f=radon_usfftadj(R,theta,epsilon,os)
N=size(R,1);
%%adjoint for bin
if (os==2)
    tmp = zeros(os*N,size(theta,2));
    tmp(1:2:end,:) = R; 
    tmp(2:2:end,:) = R;
    R = tmp*0.5;
end
%polar grid
rho=(-N*os/2:N*os/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);
%parameters for usfft
mu=-log(epsilon)/(2*N^2);
Te=1/pi*sqrt(-mu*log(epsilon)+(mu*N)^2/4);
M=ceil(2*N*Te);
%FFT1D
F=fftshift(fft(ifftshift(R)))/N;
Fe=zeros(2*N);
for i1=0:2*M
    ell1=floor(2*N*y)-M+i1;
    for i0=0:2*M
        ell0=floor(2*N*x)-M+i0;
        w0=ell0/(2*N)-x;
        w1=ell1/(2*N)-y;
        w=pi/mu*exp(-pi*pi/mu*(w0.*w0+w1.*w1));                
        %wrap
        ind1 = mod(ell1+N+6*N,2*N);%prefer to take mod for positive numbers (otherwise different behavior for python and C++)  
        ind0 = mod(ell0+N+6*N,2*N) ;               
        Fe = Fe+accumarray([ind1(:)+1 ind0(:)+1],F(:).*w(:),[2*N 2*N]);
    end
end
%IFFT2D
fe=fftshift(ifft2(fftshift(Fe)));
%unpadding and division by phi
[xx,yy]=ndgrid((-N/2:N/2-1),(-N/2:N/2-1));
phi=exp(-mu*(xx.*xx+yy.*yy));
f=fe(end/2+1+(-N/2:N/2-1),end/2+1+(-N/2:N/2-1))./phi;

