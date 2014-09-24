function [f,b_bestparam,a_bestparam gbests] = pso(x, y)
%function [f,b_bestparam,a_bestparam, traj] = pso(x, y)
%Particle Swarm Optimization algorithm for determining the filter coefficients in the filtering network between x and y
%
% Outputs:
% f: estimate for y
% b, a: filter coefficients (see filter)
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

%TODO: bdim, adim are hardcoded
bdim = 4;
adim = 3;
zdim = size(x,1);
numparticles = 20;

%pso constants
c1 = 2;
c2 = 2;
w = 0.8;
vpartlimit = 0.2;

%initializing best values;

%pbest values
pbest = ones(numparticles,1).*inf;
%pbest parameters
b_pbest = zeros(numparticles,bdim,zdim);
a_pbest = zeros(numparticles,adim,zdim);
%gbest parameters
b_gbest = zeros(numparticles,bdim,zdim);
a_gbest = zeros(numparticles,adim,zdim);

%Choosing the first numparticles from multiplier*numparticles
multiplier = 50;

c = [];
%evaluate all as and bs
values = 0;
newas = [];
newbs = [];
while (values<numparticles)
	bs = 2*rand(multiplier*numparticles,bdim,zdim)-1;
	as = 2*rand(multiplier*numparticles,adim,zdim)-1;

	max_as = size(as,1);

	for i = 1:numparticles
		if arestable(squeeze(as(i,:,:)));
			newas = [newas; as(i,:,:)];
			newbs = [newbs; bs(i,:,:)];
			values = values +1;
		end
		if (values>=numparticles) %should be ==
		    break;
		end;
	end
end
%toc;

%Replacing the actual generation with the new one

as = newas;
bs = newbs;

%random velocities
%v_bs = randn(size(bs));
%v_as = randn(size(as));

v_bs = zeros(size(bs));
v_as = zeros(size(as));

%for the best fit parameters
thebest = inf;
b_bestparam = zeros(bdim,zdim);
a_bestparam = zeros(adim,zdim);

numiters = 1000;

gbests = inf(1, numiters);
gbest = inf;
for iter = 1:numiters
for i = 1:numparticles

	%uncomment this to tract the trajectories
    %if i == 1
%	traj(iter,:) = as(i,:,1);
%    end

	%traj(:,:,iter) = particles;

	%evaluating actual particles
	f = multifilter(squeeze(bs(i,:,:)),squeeze(as(i,:,:)),x);

	%ccf = corrcoef(y,f);
	partfit = mse(y-f);%*(1-ccf(2,1));% sum((y-f).^2);

	%local param fit
	if partfit < pbest(i) && arestable(squeeze(as(i,:,:)))
	    pbest(i) = partfit;
	    b_pbest(i,:,:) = bs(i,:,:);
	    a_pbest(i,:,:) = as(i,:,:);

        %global param fit
        if partfit < gbest
            gbest = partfit;
            b_gbest = repmat(bs(i,:,:),numparticles,1);
            a_gbest = repmat(as(i,:,:),numparticles,1);
            if partfit < thebest
                thebest = partfit;
                b_bestparam = squeeze(bs(i,:,:));
                a_bestparam = squeeze(as(i,:,:));
            end;
        end;
    end;

	end;

    %condition to finish
    %if pbest < bestlimit
%	break;
%    end

    %update particles speed
    v_bs = w*v_bs + c1 * rand(numparticles,bdim,zdim) .* (b_pbest - bs) + c2 * rand(numparticles,bdim,zdim) .* (b_gbest - bs);
    v_as = w*v_as + c1 * rand(numparticles,adim,zdim) .* (a_pbest - as) + c2 * rand(numparticles,adim,zdim) .* (a_gbest - as);

    %total speed of particles
    vsums = sqrt(sum(sum(v_as.^2,3),2)+sum(sum(v_bs.^2,3),2));
    indexes = find(vsums>vpartlimit);
    for i = 1:length(indexes)
        v_as(indexes(i),:,:) = v_as(indexes(i),:,:) * (vpartlimit / vsums(indexes(i)));
        v_bs(indexes(i),:,:) = v_bs(indexes(i),:,:) * (vpartlimit / vsums(indexes(i)));
    end
    bs = bs + v_bs;
    as = as + v_as;

    gbests(iter) = gbest;
end;
f = multifilter(b_bestparam,a_bestparam,x);
