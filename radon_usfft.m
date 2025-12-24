function R=radon_usfft(f,theta,epsilon,os)
N=size(f,1);
%polar grid
rho=(-N*os/2:N*os/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);

%parameters for usfft
mu=-log(epsilon)/(2*N^2);
Te=1/pi*sqrt(-mu*log(epsilon)+(mu*N)^2/4);
M=ceil(2*N*Te);
%padding and division by phi
[xx,yy]=ndgrid((-N/2:N/2-1),(-N/2:N/2-1));
phi=exp(-mu*(xx.*xx+yy.*yy));
fe=zeros(2*N);
fe(end/2+1+(-N/2:N/2-1),end/2+1+(-N/2:N/2-1))=f./phi;
%FFT2D
Fe=fftshift(fft2(fftshift(fe)))/(4*N*N);
%gather
F=zeros(size(x));
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
        F=F+w.*Fe(ind1+1+2*N*ind0);
    end
end
%IFFT1D
R=fftshift(ifft(ifftshift(F)));

% bin
if (os==2)
    R = R(1:os:end,:)+R(2:os:end,:);
end

