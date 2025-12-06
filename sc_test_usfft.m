N=128;Ntheta=90;
%f=phantom(N);
f=rand(N,N)+1j*rand(N,N);
theta=(0:Ntheta-1)/Ntheta*pi;
epsilon=1e-12;%usfft accuracy
g=radon_usfft(f,theta,epsilon);
g_test=radon_usft(f,theta);

norm(g-g_test)

