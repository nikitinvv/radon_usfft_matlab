function f=radon_usfftadj(R,theta,epsilon,filter)
N=size(R,1);
%polar grid
rho=(-N/2:N/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);
%parameters for usfft
mu=-log(epsilon)/(2*N^2);
Te=1/pi*sqrt(-mu*log(epsilon)+(mu*N)^2/4);
M=ceil(2*N*Te);
%FFT1D
F=fftshift(fft(ifftshift(R)))/N;
if(filter==1)%paganin filter
    F=F.*abs(rho).*(.54 + .46 * cos(pi*rho./0.5))*pi/numel(theta);%step in angles
end
%scatter
x(find(x>=0.5))=0.5-1e-8;%border control
y(find(y>=0.5))=0.5-1e-8;
Fee=zeros(2*N+2*M,2*N+2*M);
for i1=0:2*M
    ell1=floor(2*N*y)-M+i1;
    for i0=0:2*M
        ell0=floor(2*N*x)-M+i0;
        w0=ell0/(2*N)-x;
        w1=ell1/(2*N)-y;
        w=pi/mu*exp(-pi*pi/mu*(w0.*w0+w1.*w1));                
        Fee = Fee+accumarray([N+M+ell1(:)+1 N+M+ell0(:)+1],F(:).*w(:),[(2*N+2*M) (2*N+2*M)]);
    end
end
%wrap
Fee(2*N+1:2*N+M,:)=Fee(2*N+1:2*N+M,:)+Fee(1:M,:);
Fee(M+1:2*M,:)=Fee(M+1:2*M,:)+Fee(end-M+1:end,:);
Fee(M+1:2*N+M,2*N+1:2*N+M)=Fee(M+1:2*N+M,2*N+1:2*N+M)+Fee(M+1:2*N+M,1:M);
Fee(M+1:2*N+M,M+1:2*M)=Fee(M+1:2*N+M,M+1:2*M)+Fee(M+1:2*N+M,end-M+1:end);
Fe=Fee(M+1:end-M,M+1:end-M);
%IFFT2D
fe=fftshift(ifft2(fftshift(Fe)));
%unpadding and division by phi
[xx,yy]=ndgrid((-N/2:N/2-1),(-N/2:N/2-1));
phi=exp(-mu*(xx.*xx+yy.*yy));
f=fe(end/2+1+(-N/2:N/2-1),end/2+1+(-N/2:N/2-1))./phi;
%circ support
circ=(2/N*xx).^2+(2/N*yy).^2<1-1/N;
f=f.*circ;


