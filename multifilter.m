function f = multifilter(bs, as, x)
%function f = multifilter(bs, as, x)
%relizes a filter network and f is the filtered input (x) using the coefficients bs and as
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

zdim = size(bs,2);
for j = 1:zdim
	fs(j,:) = filter(bs(:,j),[1; as(:,j)],x(j,:));
end;

f = mean(fs,1);
