function R=radon_usfft(f,theta,epsilon)
N=size(f,1);
%circ support
[xx,yy]=ndgrid((-N/2:N/2-1),(-N/2:N/2-1));
circ=(2/N*xx).^2+(2/N*yy).^2<1-1/N;
f=f.*circ;
%polar grid
rho=(-N/2:N/2-1)'/N;
x=rho*cos(theta);
y=rho*sin(theta);
x(find(x>=0.5))=0.5-epsilon;%border control
y(find(y>=0.5))=0.5-epsilon;


%parameters for usfft
mu=-log(epsilon)/(2*N^2);
Te=1/pi*sqrt(-mu*log(epsilon)+(mu*N)^2/4);
M=ceil(2*N*Te);
%padding and division by phi
phi=exp(-mu*(xx.*xx+yy.*yy));
fe=zeros(2*N);
fe(end/2+1+(-N/2:N/2-1),end/2+1+(-N/2:N/2-1))=f./phi;
%FFT2D
Fe=fftshift(fft2(fftshift(fe)))/(4*N*N);
%wrap
Fee=zeros(2*N+2*M,2*N+2*M);
Fee(M+1:end-M,M+1:end-M)=Fe;

Fee(1:M,:)=Fee(2*N+1:2*N+M,:);
Fee(end-M+1:end,:)=Fee(M+1:2*M,:);
Fee(:,1:M)=Fee(:,2*N+1:2*N+M);
Fee(:,end-M+1:end)=Fee(:,M+1:2*M);
%gather
F=zeros(size(x));
for i1=0:2*M
    ell1=floor(2*N*y)-M+i1;
    for i0=0:2*M
        ell0=floor(2*N*x)-M+i0;
        w0=ell0/(2*N)-x;
        w1=ell1/(2*N)-y;
        w=pi/mu*exp(-pi*pi/mu*(w0.*w0+w1.*w1));                
        F=F+w.*Fee(N+M+ell1+1+(2*N+2*M)*(N+M+ell0));
    end
end


%xeq = -N:N-1
%Fs=zeros(size(x));
%for k=1:size(x):
    %for i0=1:N:
     %   for i1=1:N:
            %Fs[k] = Fs[k]+f[i0, i1] * \
                %np.exp(-2*np.pi*1j*(xeq[i0]*y[k]+xeq[i1]*x[k]))

%norm(F-Fs,'fro')
%norm(F,'fro')

%IFFT1D
R=fftshift(ifft(ifftshift(F)));
