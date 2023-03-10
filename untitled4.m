

xinput = randn(1e4,1);

IRlength = 100;

IRtrue = zeros(IRlength,1);
IRtrue(3) = 1;
IRtrue(12) = -0.2;
IRtrue(19) = 0.12;
IRtrue(30) = 1;
IRtrue(43) = 1;
IRtrue(57) = -0.6;
IRtrue(79) = 0.3;
IRtrue(92) = -0.1;

youtput = contwo(xinput,IRtrue);
youtput = youtput(1:length(xinput));

youtput = youtput + 0.01*randn(size(youtput));

convfactor = 0.01;
nsteps = length(xinput);

[IRestimate,errorhistory] = LMS_Ncoeffs(xinput,youtput,IRlength,convfactor,nsteps);

figure(1)
semilogy(abs(errorhistory))
grid

figure(2)
ivec = [1:IRlength];
plot(ivec,IRtrue,ivec,IRestimate,'*');
grid






