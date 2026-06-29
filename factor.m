function[b,se]=factor(s)
%  function[b,se]=factor(s)
%
%  This program accepts input in form of spectral density s
%  of a univariate process, and returns the kernel of the Wold representation
%  for the process in b.  The one step ahead prediction
%  error variance is in se.
%  The method is described in Whittle's Prediction and Regulation
%  (2nd ed, U of Minn Press, 1983) on page 26.
%
%  We take log(s), call it ls.  Then we take the inverse fft
%  of ls, call it cc.  The one-step ahead prediction error is given
%  by the Wiener-Kolmogorov formula se=exp(cc(1)).  We take the
%  strictly positive powers of cc, pad with zeros, take the fft of
%  the result, take its exponential, then finally inverse fft to get
%  the moving average kernel, which is stored in b.  The interesting
%  outputs are se, the one-step ahead prediction variance, and
%  b, the kernel in the moving average representation.
%
%  SEE ALSO: FACTORF.M, FORECASF.M
ls=log(s);
cc=ifft(ls);
n0=max(size(cc));
n02=n0/2;
c1=cc(2:n0/2+1);
se=exp(cc(1));
p1=[0;c1;zeros(n02,1)];
pf=fft(p1);
bf=exp(pf);
bb=ifft(bf);
b=real(bb);
