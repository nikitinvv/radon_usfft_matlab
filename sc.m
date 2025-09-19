N=256;Ntheta=180;
f=phantom(N);
theta=(0:Ntheta-1)/Ntheta*pi;
epsilon=1e-12;%usfft accuracy
g=radon_usfft(f,theta,epsilon);
ff=radon_usfftadj(g,theta,epsilon,1);
imagesc(real([f ff]));

%adjoint test
fwd=@(f)radon_usfft(f,theta,epsilon);
adj=@(g)radon_usfftadj(g,theta,epsilon,0);
f=rand(N);
g=rand(N,Ntheta);
gg=fwd(f);
ff=adj(g);

(sum(conj(ff(:)).*f(:))-sum(gg(:).*conj(g(:))))/sum(conj(ff(:)).*f(:))
