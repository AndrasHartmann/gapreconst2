function [ff q1 q2 gbests] = reconstruct(st, slength)
% function [ff q1 q2] = reconstruct(st, slength)
% Reconstruct the record with identifier st where fitting is made of slength datapoints on the prior
%
% Return:
%   ff: Reconstructed time-series
%   q1: Q1 error as defined in [1]
%   q2: Q2 error as defined in [1]
% Connected publications:
% 
% suppose tests are downloaded into data folder
% The dataset can be downloaded by using the wfdb package from Physionet
% $wfdb2mat -r challenge/2010/set-a/a11
%
%    GapReconst2, reconstruction of physiological signals
%    Copyright (C) 2014 Andras Hartmann <hdbandi@gmail.com>
%
%    If you use this program for scientific publishing,
%    please cite the following paper:
%
%    Hartmann, A., Lemos, J. M., Costa, R. S., & Vinga, S. (2014). 
%    Identifying IIR filter coefficients using particle swarm optimization 
%    with application to reconstruction of missing cardiovascular signals.
%    Engineering Applications of Artificial Intelligence, 34, 193–198.
%    doi:10.1016/j.engappai.2014.05.014
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, version 3 of the License.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%load dataset
load(['data/' st 'm.mat']);
%find what is missing
for i = 1:size(val,1)
	if val(i,end-3749:end) == zeros(1,3750)
		break;
	end;
end;
ind = i;

%determine signal names from header
myfile = fopen(['data/' st 'm.hea'], 'r');
tline = fgetl(myfile);
signalnames = {};
while 1
    tline = fgetl(myfile);
    if ~ischar(tline),
	break;
    end
    A = strread(tline,'%s');
    A = A(end);
    signalnames{end+1} = A{:};
end
fclose(myfile);

%c = corrcoef(val');
%el = find(abs(c(ind,:)) > 0.1 & abs(c(ind,:)<1));
x = val(:,end-slength:end-3750);
xx = val(:,end-3749:end);
x(ind,:) = [];
xx(ind,:) = [];

if strcmp(signalnames{ind},'RESP')
    %trick: filter for low freqz (use only for resp)
    load filtcoeff;

    %freqz(bl,al,1024,Fs)
    x = filtfilt(bl,al,x')';
    xx = filtfilt(bl,al,xx')';
end

%remove mean and std
mx = mean(x,2);
mxx = mean(xx,2);
x = x - repmat(mx,1,size(x,2));
%x = x ./ repmat(sx,1,size(x,2));
xx = xx - repmat(mxx,1,size(xx,2));
%xx = xx ./ repmat(sxx,1,size(xx,2));
y = val(ind,end-slength:end-3750);
my = mean(y);
sy = std(y);

%check if the function is not identical
if (sy == 0)
    ff = ones(1,3750)*my;
    q1 =1; q2=1; 
    return;
end
%TODO: we should rather exclude identical signals
%well, this is tricky, but still :)

% subtract the mean
y = y-my;
%y = y./sy;
%run genetic algo
%[f,b,a] = genalg(x, y);
[f,b,a,gbests] = pso(x, y);
ff = multifilter(b, a, xx);
yy = load(['data/' st '.missing']);
%filtval
%ff = ff*sy;
%add the mean again
ff = ff+my;

ccf = corrcoef(yy,ff);
q1 = 1-mse(yy-ff')/var(yy);
q2 = ccf(2);

