clc;
clear all;
close all;

xs = dlmread('xs.txt');
xs_in = dlmread('xs_in.txt');
ys = dlmread('ys.txt');
ys_in = dlmread('ys_in.txt');
alpha = dlmread('alpha_s.txt');
alpha_in = dlmread('alpha_s_in.txt');
b = dlmread('b.txt');
train_class1_data = dlmread('train_class1_data.txt');
train_class2_data = dlmread('train_class2_data.txt');
test_class1_data = dlmread('test_class1_data.txt');
test_class2_data = dlmread('test_class2_data.txt');

% h=0.2;
% points = -1:h:1;
% N = length(points);
% label = ones(N^4,1);
% p=1;
% for i=1:N
%     for j=1:N
%         for k=1:N
%             for l=1:N
%                 feature_vec = [points(i), points(j), points(k), points(l)];
%                 label(p) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%                 p = p+1;
%             end
%         end
%     end
% end

h=0.2;
points = -1:h:1;
N = length(points);

figure(1);
% subplot(231);
test_label = ones(N,N);
true_label = ones(N,N);
for i=1:N
    for j=1:N
        feature_vec = [points(i), points(j)];
        true_label(i,j) = label_func(feature_vec);
        test_label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, b, feature_vec);
    end
end
[X, Y] = meshgrid(points,points);
contourf(X,Y,true_label); %x,y (z=0, w=0)
err = abs(true_label-test_label);
[c,temp] = size(err(err==0));
accuracy = c*100/N/N
% contourf(X,Y,test_label); %x,y (z=0, w=0)
% hold on;
% plot(X,Y,'ro');

% 
% subplot(232);
% label = ones(N,N);
% for i=1:N
%     for j=1:N
%         feature_vec = [points(i), 0, points(j), 0];
%         label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%     end
% end
% [X, Y] = meshgrid(points,points);
% contourf(X,Y,label); %x,z (y=0, w=0)
% 
% subplot(233);
% label = ones(N,N);
% for i=1:N
%     for j=1:N
%         feature_vec = [points(i), 0, 0, points(j)];
%         label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%     end
% end
% [X, Y] = meshgrid(points,points);
% contourf(X,Y,label); %x,w (y=0, z=0)
% 
% subplot(234);
% label = ones(N,N);
% for i=1:N
%     for j=1:N
%         feature_vec = [0, points(i), points(j), 0];
%         label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%     end
% end
% [X, Y] = meshgrid(points,points);
% contourf(X,Y,label); %y,z (x=0, w=0)
% 
% subplot(235);
% label = ones(N,N);
% for i=1:N
%     for j=1:N
%         feature_vec = [0, points(i), 0, points(j)];
%         label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%     end
% end
% [X, Y] = meshgrid(points,points);
% contourf(X,Y,label); %y,w
% 
% subplot(236);
% label = ones(N,N);
% for i=1:N
%     for j=1:N
%         feature_vec = [0, 0, points(i), points(j)];
%         label(i,j) = test(xs, xs_in, ys, ys_in, alpha, alpha_in, feature_vec);
%     end
% end
% [X, Y] = meshgrid(points,points);
% contourf(X,Y,label); %z,w
