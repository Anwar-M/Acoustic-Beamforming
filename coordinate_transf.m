z_offset = -0.1;

P = [-1628 -4054 -8061;
      2370 -4034 -8064;
      2386 -6558 -6438;
     -1611 -6598 -6469]';
 
X = [0 0 1; 4 0 1; 4 3 1; 0 3 1]';


%% Check, do I get my original coordinates back?

C = P(:,2:4)/X(:,2:4);

C*X(:,2:4)

%% We have fourth point, how large is the error when comparing this point
%  compared with calculated value?

z = 1:4;

for j = 1:4
    A = P(:,z(z~=j))/X(:,z(z~=j));
    (A*X(:,j)-P(:,j))
end

%%
P = P';
D = det(P(1:3,:));
d = 1;

a = -d*det([ones(3,1) P(1:3,2) P(1:3,3)])/D;
b = -d*det([P(1:3,1) ones(3,1) P(1:3,3)])/D;
c = -d*det([P(1:3,1) P(1:3,2) ones(3,1)])/D;

%% Add offset of 10 cm in z direction local coordinates

a*P(4,1)+b*P(4,2)+c*P(4,3) + d;

[Xs,Ys] = meshgrid(-1620:500:2380, -6600:500:-4040);
Zs = (-d-a*Xs-b*Ys)/c;
figure; 
hold on;
surf(Xs,Ys,Zs);

norm_vect = [a; b; c];
norm_vect = norm_vect/norm(norm_vect);

Xs = Xs + 100*norm_vect(1);
Ys = Ys + 100*norm_vect(2);
Zs = Zs + 100*norm_vect(3);

surf(Xs,Ys,Zs);hold off;

P;
PP = P + [[100*ones(4,1)*norm_vect(1)] [100*ones(4,1)*norm_vect(2)] [100*ones(4,1)*norm_vect(3)]];


%% Transform all microphone positions to tunnel coordinates