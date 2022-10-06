function [h_estimate,errorhistory] = MC_LMS_Ncoeffs(xinput,youtput,n_sources,n_measurements,h_length,convfactor,n_steps,h_startestimate)
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
%   n_sources   Number of sources in xinput
%   n_measurements  Number of measurements in youtput
%   h_length    The desired length of the h_estimate, in samples.
%   convfactor  The convergence factor. Try a value like 0.1 or smaller.
%   nsteps
%   h_startestimate  (optional) A start h estimate
%
% Output parameters:
%   h_estimate  The LMS estimate of h, of size [h_length,1]
%   errorhistory    A vector, of size [1,ninput], containing the error
%                   after each iteration step
%
% peter.svensson@ntnu.no 2 June 2022
%
% [h_estimate,errorhistory] = LMS_Ncoeffs(xinput,youtput,h_length,convfactor,nsteps,h_startestimate);

if nargin < 8
    h_estimate = zeros(h_length,1);
else
    if length(h_startestimate) ~= h_length
        error('ERROR: An IRstartestimate must have the length = IRlength')
    end
    h_estimate = h_startestimate;
end

for i = 1:n_measurements
    alpha = 1/mean(xinput(:,i).^2)*convfactor;

    nsig = length(xinput);
    errorhistory = zeros(1,nsig);

    for ii =1:n_steps
        flippedx = xinput(ii:-1:max([1 ii-h_length+1]),i);
        lengthflippedx = length(flippedx);
        if lengthflippedx < h_length
            flippedx = [flippedx;zeros(h_length-lengthflippedx,1)];
        end

        y_hat = sum(flippedx.*h_estimate);

        e =  youtput(ii) - y_hat ;
        errorhistory(ii) = e;

        h_estimate = h_estimate + alpha*e*flippedx;
    end
end
