function verify_code
% Independent derivative, dense-matrix, and conditional-GCV checks.
folder=fileparts(mfilename('fullpath')); previous_path=path;
cleanup=onCleanup(@() path(previous_path)); %#ok<NASGU>
addpath(folder,'-begin');
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
report=struct;
families={'gaussian','defocus','motion'};
functions={@gauss_psf_symbol_derivs,@defocus_psf_symbol_derivs,@motion_psf_symbol_derivs};
points={[2.3;3.1],2.7,[8.3;.13]};
n=16; rng(204,'twister'); bhat=fft2(randn(n));
[u,v]=ndgrid(0:n-1); s=sqrt(4*sin(pi*u/n).^2+4*sin(pi*v/n).^2);
report.derivatives=struct('family',{},'first',{},'second',{},'gradient',{});
for k=1:3
    y=points{k}; symbol=@(z) functions{k}(z,n);
    [c,dc,d2c]=symbol(y); de=0; d2e=0; ge=0;
    [~,g]=reconstruction_objective(y,bhat,s,.12,.07,.03,ones(n),symbol,@(z) z(1)>0);
    for j=1:numel(y)
        h=1e-5*max(1,abs(y(j))); step=zeros(size(y)); step(j)=h;
        [cp,dcp]=symbol(y+step); [cm,dcm]=symbol(y-step);
        de=max(de,norm((cp-cm)/(2*h)-dc(:,j))/max(1,norm(dc(:,j))));
        d2e=max(d2e,norm((dcp-dcm)/(2*h)-d2c(:,:,j),'fro')/max(1,norm(d2c(:,:,j),'fro')));
        fp=reconstruction_objective(y+step,bhat,s,.12,.07,.03,ones(n),symbol,@(z) z(1)>0);
        fm=reconstruction_objective(y-step,bhat,s,.12,.07,.03,ones(n),symbol,@(z) z(1)>0);
        ge=max(ge,abs((fp-fm)/(2*h)-g(j))/max(1,abs(g(j))));
    end
    assert(max([de,d2e,ge])<1e-5 && abs(c(1)-1)<1e-12);
    report.derivatives(k)=struct('family',families{k},'first',de,'second',d2e,'gradient',ge);
    fprintf('%s derivative errors: %.3g %.3g %.3g\n',families{k},de,d2e,ge);
end

% Independent small dense matrices check FFT normalization and the GCV trace.
n=8; N=n^2; eta=37; delta=2; lambda=sqrt(delta/eta);
[u,v]=ndgrid(0:n-1); s=sqrt(4*sin(pi*u/n).^2+4*sin(pi*v/n).^2);
symbol=@(z) gauss_psf_symbol_derivs(z,n); y=[1.2;.8];
c=reshape(symbol(y),n,n); A=zeros(N); LTL=zeros(N);
for j=1:N
    e=zeros(n); e(j)=1;
    col=real(ifft2(c.*fft2(e))); A(:,j)=col(:);
    col=real(ifft2(s.^2.*fft2(e))); LTL(:,j)=col(:);
end
b=reshape(sin((1:N)/7),n,n); bhat=fft2(b);
Q=eta*(A'*A)+delta*LTL; mu=Q\(eta*A'*b(:)); R=chol(Q);
marginal=eta/2*norm(A*mu-b(:))^2+delta/2*(mu'*LTL*mu)+sum(log(diag(R)));
f=reconstruction_objective(y,bhat,s,lambda,1/sqrt(eta),delta/eta,s,symbol,@(z) all(z>0));
report.marginal_identity_error=abs(marginal-eta*f-N/2*log(eta));
hat=A*((A'*A+lambda^2*LTL)\A');
matrix_gcv=(norm(b(:)-hat*b(:))^2/N)/(1-trace(hat)/N)^2;
report.gcv_error=abs(matrix_gcv-conditional_gcv(lambda,c,s,bhat));
assert(report.marginal_identity_error<1e-9 && report.gcv_error<1e-10);

fprintf('PASS: PSF derivatives, objective gradients, dense marginalization, and GCV.\n');
verify_solver_inputs;
end
