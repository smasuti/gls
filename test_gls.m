% This script tests the general least squares. 
% Author : Sagar Masuti
% Date   : 04-02-2018
% -------------------------------------------------------------------------

%% First Example : Pearson's data
filename = 'pearson01.txt';
d=load (filename);

param=ones(1,length(d(:,2)));
G=[d(:,2),param'];
m=inv(G'*G)*G'*d(:,4);

[A,S]=gls(m', d(:,2),d(:,4),[d(:,3) d(:,5)]);

figure(2);clf;
errorbarxy(d(:,2),d(:,4),d(:,3),d(:,5),{'ko', 'b', 'r'});
hold on;
m=A(:,1)*d(:,2)+A(:,2);
fprintf ('m  : %.2e +/- %.2e\n',A(:,1),sqrt(S(1,1)));
fprintf ('b  : %.2e +/- %.2e\n',A(:,2),sqrt(S(2,2)));
plot(d(:,2),m,'k','LineWidth',1.5);
axis tight
xlabel('X');
ylabel('Y');

%% Second Example: 1 dependent variable and 1 independent variables
x_err=[0.8636,0.9430,0.2834,0.0421,0.0748,0.2833,0.0566,0.2079,0.4753,0.5260];
y_err=[0.8636,0.9430,0.7381,0.5446,0.0384,0.1331,0.0793,0.8751,0.8462,0.4990];
x = (1:10); 
y = 8*(x')+ 10 + y_err' ; 
[A,S]=anisotropy.gls([7 8], x',y,[x_err' y_err']);

figure(3);clf;
errorbarxy(x',y,x_err',y_err',{'ko', 'b', 'r'});
hold on;
m=A(:,1)*x'+A(:,2);
fprintf ('m  : %.2e +/- %.2e\n',A(:,1),(S(1,1)));
fprintf ('b  : %.2e +/- %.2e\n',A(:,2),(S(2,2)));
plot(x,m,'k','LineWidth',1.5);
axis tight
xlabel('X');
ylabel('Y');

%% Third Example :  1 dependent variable and 2 independent variables.
x = [0.8425,0.7011,0.0659,0.2509,0.0353,0.0837,0.8020,0.0591,0.6363,0.7515];
z = [0.7317,0.0692,0.9226,0.0438,0.4857,0.8959,0.6124,0.3552,0.4364,0.8492];
y_err=[-0.0880,-0.0278,0.1546,0.0419,-0.0656,0.1014,0.0587,0.0159,-0.1804,-0.2086];

% Synthetic data
y = (3-2*x-5*z)+y_err; % Equation of the plane : 2*x+y+5*z=3

param=ones(1,length(x'));
G=[x',z',param'];
w=diag((1./y_err).^2);
m=inv(G'*w*G)*G'*w*y';
cm=var(y)*inv(G'*w*G);
res=sqrt(sum((m(1)*x+m(2)*z+m(3)-y).^2));

fprintf ('-------- Using weighted least squares -----------\n')
fprintf ('residual  : %.2e\n',res);
fprintf ('a  : %.2e +/- %.2e\n',m(1),sqrt(cm(1,1)));
fprintf ('b  : %.2e +/- %.2e\n',m(2),sqrt(cm(2,2)));
fprintf ('c  : %.2e +/- %.2e\n',m(3),sqrt(cm(3,3)));

x_err=0.1*[0.9510,-1.0733,0.1884,1.5951,-0.7156,-0.8268,-0.8024,-0.1677,0.4377,0.6016];
z_err=0.1*[-0.5256,-1.3830,-1.9557,-0.3622,-0.6027,1.3365,0.1444,-0.6664,-0.8471,0.4049];

fprintf ('-------- Using General least squares -----------\n')
[A,S]=gls(m', [x' z'],y',[x_err' z_err' y_err']);
fprintf ('a  : %.2e +/- %.2e\n',A(:,1),sqrt(S(1,1)));
fprintf ('b  : %.2e +/- %.2e\n',A(:,2),sqrt(S(2,2)));
fprintf ('c  : %.2e +/- %.2e\n',A(:,3),sqrt(S(3,3)));

