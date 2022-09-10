
xinput = randn(1e2,1);

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

[xinputestimate,errorhistory] = LMS_Ncoeffs(IRtrue,youtput,1e2,convfactor,nsteps);
totalerrorhistory = errorhistory;
i = 1000;

while i > 0
    [xinputestimate,errorhistory] = LMS_Ncoeffs(IRtrue,youtput,1e2,convfactor,nsteps,xinputestimate);
    totalerrorhistory = [totalerrorhistory errorhistory];
    i = i-1;
end
figure(1)
semilogy(abs(errorhistory))
grid

figure(2)
ivec = [1:IRlength];
plot(ivec,xinput,ivec,xinputestimate,'*');
grid

%[xinputestimate,errorhistory] = LMS_Ncoeffs(IRtrue,youtput,1e2,convfactor,nsteps, xinputestimate);






