function H = computeHomography(pts1, pts2)
% Compute homography that maps from pts1 to pts2 using least squares solver
% 
% Input: pts1 and pts2 are 3xN matrices for N points in homogeneous
% coordinates. 
%
% Output: H is a 3x3 matrix, such that pts2~=H*pts1


% Normalized
pts1_transpose = transpose(pts1);
pts2_transpose = transpose(pts2);

S1 = std(pts1_transpose);
S2 = std(pts2_transpose);

T1 = zeros(3, 3);
T1(1, 1) = 1/S1(1);
T1(2, 2) = 1/S1(2);
T1(3, 3) = 1;
pts1_normalized = pts1_transpose*T1;
pts1_normalized = transpose(pts1_normalized);

T2 = zeros(3, 3);
T2(1, 1) = 1/S2(1);
T2(2, 2) = 1/S2(2);
T2(3, 3) = 1;
pts2_normalized = pts2_transpose*T2;
pts2_normalized = transpose(pts2_normalized);

% Compute homography by DLT

% form matrix A
A = zeros(4, 9);

for i = 1 : 4
    A(i*2-1, 1) = -pts1_normalized(1, i);
    A(i*2-1, 2) = -pts1_normalized(2, i);
    A(i*2-1, 3) = -1;
    
    A(i*2-1, 7) = pts1_normalized(1, i)*pts2_normalized(1, i);
    A(i*2-1, 8) = pts1_normalized(2, i)*pts2_normalized(1, i);
    A(i*2-1, 9) = pts2_normalized(1, i);
    
    A(i*2, 4) = -pts1_normalized(1, i);
    A(i*2, 5) = -pts1_normalized(2, i);
    A(i*2, 6) = -1;
    
    A(i*2, 7) = pts1_normalized(1, i)*pts2_normalized(2, i);
    A(i*2, 8) = pts1_normalized(2, i)*pts2_normalized(2, i);
    A(i*2, 9) = pts2_normalized(2, i);    
end

[U, S, V] = svd(A);
h = V(:, end);
% convert h to 3x3 matrix H
H = zeros(3, 3);
for i = 1:3
    H(i, 1) = h(i*3-2);
    H(i, 2) = h(i*3-1);
    H(i, 3) = h(i*3);
end

% Unnormalize H
H = inv(T2)*H*T1;



