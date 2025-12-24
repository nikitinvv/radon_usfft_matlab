N=256;Ntheta=N*2;
f = zeros(N)+1j*0;
f(N/8+1:N/8+N-N/4,N/8+1:N/8+N-N/4)=phantom(N-N/4);

theta=(0:Ntheta-1)/Ntheta*pi;
epsilon=1e-3;%usfft accuracy
os = 2;% oversampling in frequencies

% generate data
d=radon_usfft(f,theta,epsilon,os);


minf=@(Ru,d)norm(Ru - d)^2;
R=@(f)radon_usfft(f,theta,epsilon,2);
RT=@(g)radon_usfftadj(g,theta,epsilon,2);
redot=@(f,g)real(dot(f(:), g(:)));
u = f*0;

Ru = R(u);
disp('iterations')
for k=0:128
    grad = RT(2 * (Ru - d));
    Rgrad = R(grad);
    if k == 0
        eta = -grad;
        Reta = -Rgrad;
    else
        beta = redot(Rgrad, Reta) / redot(Reta, Reta);
        eta = beta * eta - grad;
        Reta = beta * Reta - Rgrad;
    end
    alpha = -redot(grad, eta) / (2 * redot(Reta, Reta));
    u = u+alpha * eta;
    Ru = Ru+alpha * Reta;
    if mod(k,8) == 0
        fprintf("%d, %f\n", k, minf(Ru, d));
    end
end

%figure(1);imagesc([real(u) imag(u)]);colorbar;
figure(2);imagesc([real(u)-real(f)]);colorbar;
      