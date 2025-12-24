N=256;Ntheta=N/4;

theta=(0:Ntheta-1)/Ntheta*pi;
epsilon=1e-12;%usfft accuracy

g=radon_usfft(f,theta,epsilon,1); %%% os = 1
g2=radon_usft(f,theta,1); %%% os = 1
norm(g-g2)

g=radon_usfft(f,theta,epsilon,2); %%% os = 2
g2=radon_usft(f,theta,2); %%% os = 2
norm(g-g2)

%adjoint test
fwd=@(f)radon_usfft(f,theta,epsilon,1);
adj=@(g)radon_usfftadj(g,theta,epsilon,1);
f=rand(N);
g=rand(N,Ntheta);
gg=fwd(f);
ff=adj(g);

(sum(conj(ff(:)).*f(:))-sum(gg(:).*conj(g(:))))/sum(conj(ff(:)).*f(:))

%adjoint test
fwd=@(f)radon_usfft(f,theta,epsilon,2);
adj=@(g)radon_usfftadj(g,theta,epsilon,2);
f=rand(N);
g=rand(N,Ntheta);
gg=fwd(f);
ff=adj(g);

(sum(conj(ff(:)).*f(:))-sum(gg(:).*conj(g(:))))/sum(conj(ff(:)).*f(:))
