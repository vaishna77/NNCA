clc;
clear all;
close all;

Qchoice = input('Enter 0 for Matern and 1 for Gaussian: ');
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
N1 = length(test_class1_data);

test_data_label = dlmread('test_label.txt');
test_data_label1 = test_data_label(N1+1:length(test_data_label),:);
test_data_label2 = test_data_label(1:N1,:);

h=0.01;
if (Qchoice == 0)
    points_x = -1.0:h:1.0;
else
    points_x = -1.4:h:1.4;
end
N_x = length(points_x);
points_y = -1.0:h:1.0;
N_y = length(points_y);

figure(1);
test_label = ones(N_y,N_x);
true_label = ones(N_y,N_x);
for i=1:N_x
    for j=1:N_y
        feature_vec = [points_x(i), points_y(j)];
        true_label(j,i) = label_func(feature_vec);
        [out, test_label(j,i)] = test(xs, xs_in, ys, ys_in, alpha, alpha_in, b, feature_vec, Qchoice);
    end
end
[X, Y] = meshgrid(points_x,points_y);
contourf(X,Y,test_label);
pbaspect([1 1 1])
colormap ([0 0.4470 0.7410; 0.8500 0.3250 0.0980]);
hold on
% err = abs(true_label-test_label);
% [c,temp] = size(err(err==0));
% accuracy = c*100/N_x/N_y
h1 = plot(train_class1_data(:,1), train_class1_data(:,2), 'ko', 'MarkerFaceColor', 'r');
h2 = plot(train_class2_data(:,1), train_class2_data(:,2), 'ko', 'MarkerFaceColor', 'b');
legend([h1 h2],{'class 1','class 2'})

figure(2);
contourf(X,Y,test_label);
pbaspect([1 1 1])
colormap ([0 0.4470 0.7410; 0.8500 0.3250 0.0980]);
hold on
h1 = plot(test_class1_data(:,1), test_class1_data(:,2), 'ko', 'MarkerFaceColor', 'r');
h2 = plot(test_class2_data(:,1), test_class2_data(:,2), 'ko', 'MarkerFaceColor', 'b');
legend([h1 h2],{'class 1','class 2'})
