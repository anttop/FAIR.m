% Copyright 2019 Lukas F. Lang and Sebastian Neumayer
%
% This file is part of TBIR.
%
%    TBIR is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    TBIR is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with TBIR.  If not, see <http://www.gnu.org/licenses/>.
function tests = addnoiseTest
    tests = functiontests(localfunctions);
end

function resultTest(testCase)

x = rand(30, 40, 50);
y = addnoise(x, 0.05);
verifyEqual(testCase, size(y), size(x));

x = randi(255, 10, 20);
y = addnoise(x, 0.05);
verifyEqual(testCase, size(y), size(x));

end
