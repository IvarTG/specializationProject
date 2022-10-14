function [IRestimate,errorhistory] = LMS_Ncoeffs(xinput,youtput,IRlength,convfactor,nsteps,IRstartestimate)
% LMS_Ncoeffs estimates an IR from an input signal and an output signal using
% the iterative LMS to find the IR estimate. The input signal must be a
% noise sequence, possibly band-pass filtered.
%
% Input parameters:
%   xinput      The input sequence, x, of size [ninput,1];
%   youtput     The output sequence, y, of size [ninput,1];
%               Note that y must be of the same length as x.
%               This means that if conv has been used to generate a longer
%               y than x, then the end of y must be cut off.
%   IRlength    The desired length of the IRestimate, in samples.
%   convfactor  The convergence factor. Try a value like 0.1 or smaller.
%   nsteps
%   IRstartestimate  (optional) A start IR estimate
%
% Output parameters:
%   IRestimate  The LMS estimate of the IR, of size [IRlength,1]
%   errorhistory    A vector, of size [1,ninput], containing the error
%                   after each iteration step
%
% peter.svensson@ntnu.no 2 June 2022
%
% [IRestimate,errorhistory] = LMS_Ncoeffs(xinput,youtput,IRlength,convfactor,nsteps,IRstartestimate);

if nargin < 6
    IRestimate = zeros(IRlength,1);
else
    if length(IRstartestimate) ~= IRlength
        error('ERROR: An IRstartestimate must have the length = IRlength')
    end
    IRestimate = IRstartestimate;
end

alpha = 1/mean(xinput.^2)*convfactor;

nsig = length(xinput);
errorhistory = zeros(1,nsig);



y_hat = contwo(xinput,IRestimate);

e =  youtput - y_hat ;
errorhistory = e;

IRestimate = IRestimate + alpha*e*(xinput.');
