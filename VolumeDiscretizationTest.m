%% VolumeDiscretizationTest
%
function tests = VolumeDiscretizationTest
    tests = functiontests(localfunctions);
end

function testIndexing(testCase)
% Verify index conversion is symmetric

% Construct volume
vd = VolumeDiscretization();
% Get dims
dims = size(vd.V);
% Generate random offset
offset = ceil(rand(size(dims)).*dims).';

% Convert offset to index and back
idx = vd.offset2index(offset);
verifyEqual(testCase,vd.index2offset(idx),offset);

% Convert offset to coordinate and back
X = vd.offset2coordinate(offset);
verifyEqual(testCase,vd.coordinate2offset(X),offset);
end