
%
% Copyright (c) 2014, Dawei Wang
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Dawei Wang, nor the names of its
%       contributors may be used to endorse or promote products derived from
%       this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

% read the images and save them in a matrix.
fid = fopen('all_train.list'); faces = [];
while !feof(fid)
  f = imread(fgetl(fid));
  faces = [faces; double(f(:)')];
end

% note that we have far more features than samples, it is more computationally feasible
% to only find out the singular value decomposition for the covariance matrix instead of the original.
[u, s, d] = svd(cov(faces)); 

% project kawamura into the new space - though we don't need all the scores but for now - 
% let's compute everything and extract some of them in the loop. then scale and plot it, 
% don't forget to add the mean and compute the reconstruction error and display it.
reduced_kawamura = center(faces(1, :)) * d; index = 0; errors = [];
for n = 1: 960
  kawamura = reduced_kawamura(:, 1: n) * d(:, 1: n)' + mean(faces(1, :));
  errors = [errors, sum(sumsq(kawamura - faces(1, :))) / prod(size(kawamura))];
end


