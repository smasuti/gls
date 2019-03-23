function [A,S]=gls(a,inx,iny,err)
% General Least squares.
% Author : Sagar Masuti
% Date   : 02-02-2018
% equation:  y = a0*x + a1*z + a2 +...
% a : vector of initial guess.

i=0;
tol=1e-16;
change=1;

while((i<1000) && (change > tol))
    newa=[a(1:end-1) 1];
    acell=repmat({newa},1,length(inx));
    B=blkdiag(acell{:});

    % inx = [x1 z1 ...;
    %        x2 z2 ...];
    J=-[inx,ones(size(inx,1),1)];

    % err = [xerr1 zerr1 ... yerr1;
    %        xerr2 zerr2 ... yerr2];
    varerr = err';
    Q=diag(varerr(:).^2);

    % iny=[y1
    %      y2];
    K = ((a*(-J'))'-iny) ;

    %            -1
    % We=(B Q B')
    We = inv(B * Q * B');

    %                          -1 
    % Solution is X = (J' We J)  J' We K

    X = inv(J'*We*J)*J'*We*K;
    Ve=J*X-K;

    S = inv(J'*We*J);
    newa=a+X';
    i=i+1;
    change=abs(min(X));
    a = newa;
    res=sum((a*[inx ones(length(inx),1)]'-iny').^2);
%     fprintf('residual : %.2e\n',res);
end
fprintf('Number of iterations: %d, change= %.2e, res=%.2e\n',i,change,res);
A=newa;
end
