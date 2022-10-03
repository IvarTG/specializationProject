
%In this code we try out using a measured IR as our input signal x, the
%systemIR h_0 is the unknown sourcesig q_0 that we want to find. Our
%estimate of the systemIR, h, is the estimate of our sourcesig, q, that we
%are updating adaptively

q_0 = randn(1e2,1);

IRlength = 100;

IR = zeros(IRlength,1);
IR(3) = 1;
IR(12) = -0.2;
IR(19) = 0.12;
IR(30) = 1;
IR(43) = 1;
IR(57) = -0.6;
IR(79) = 0.3;
IR(92) = -0.1;

youtput = contwo(q_0,IR);
youtput = youtput(1:length(q_0));

youtput = youtput + 0.01*randn(size(youtput));

convfactor = 0.01;
nsteps = length(q_0);

[q,errorhistory] = LMS_Ncoeffs(IR,youtput,1e2,convfactor,nsteps);
totalerrorhistory = errorhistory;
i = 300;

while i > 0
    [q,errorhistory] = LMS_Ncoeffs(IR,youtput,1e2,convfactor,nsteps,q);
    totalerrorhistory = [totalerrorhistory errorhistory];
    i = i-1;
end
figure(1)
semilogy(abs(errorhistory))
grid

figure(2)
ivec = [1:IRlength];
plot(ivec,q_0,ivec,q,'*');
grid

%[xinputestimate,errorhistory] = LMS_Ncoeffs(IRtrue,youtput,1e2,convfactor,nsteps, xinputestimate);






