close all
set(0,'DefaultFigureWindowStyle','docked')
lw = 1.5; fz = 24;
%
% plot convergence test - F-256
%
figure; 
load errors_F_256N;
[mf256 n] = size(errors_F_256N);
q2 = ones(mf256+1,1); q1 = q2;

N_F256 = 2*ones(mf256+1,1);
for i=1:mf256
  N_F256(i+1) = 2*N_F256(i);
  s_h = N_F256(i);
  q2(i) = 1/s_h^2*4;
  q1(i) = 1/s_h^1*2;
end

loglog(N_F256(1:mf256),errors_F_256N(:,1)./errors_F_256N(1,1),'r*--'),hold on % error
loglog(N_F256(1:mf256),errors_F_256N(:,2)./errors_F_256N(1,2),'go-.'),hold on % residual
loglog(N_F256(1:mf256),q2(1:mf256),'k--'),hold on
loglog(N_F256(1:mf256),q1(1:mf256),'b-.'),hold on
set(gca,'XTick',N_F256(1:mf256))
  legend('Error','Residual','Quadradic','Linear',3)
V=axis;
V(1)=N_F256(1);
V(2)=N_F256(mf256);
V(3)=errors_F_256N(mf256,1);
V(4)=1.5;
axis(V);

xlabel('N cells X-Y-Z direction');
ylabel('|x_N|_{\inf}/|x_0|')
title('Convergence: 1 F-cycle, constant coef. Laplacian, u=(x^4-x^2) on (0,1)^3, N=256/PE')
grid
print(gcf,'-djpeg100','converg_Fcycles_256N')
%
% plot convergence test - F-32
%
figure; 
load errors_F_032N;
[mf032 n] = size(errors_F_032N);
q2 = ones(mf032+1,1); q1 = q2;

N_F032 = 2*ones(mf032+1,1);
for i=1:mf032
  N_F032(i+1) = 2*N_F032(i);
  s_h = N_F032(i);
  q2(i) = 1/s_h^2*4;
  q1(i) = 1/s_h^1*2;
end

loglog(N_F032(1:mf032),errors_F_032N(:,1)./errors_F_032N(1,1),'r*--'),hold on
loglog(N_F032(1:mf032),errors_F_032N(:,2)./errors_F_032N(1,2),'go-.'),hold on
loglog(N_F032(1:mf032),q2(1:mf032),'k--'),hold on
loglog(N_F032(1:mf032),q1(1:mf032),'b-.'),hold on
set(gca,'XTick',N_F032(1:mf032))
legend('Error','Residual','Quadradic','Linear',1)
V=axis;
V(1)=N_F032(1);
V(2)=N_F032(mf032);
V(3)=errors_F_032N(mf032,1)/2;
V(4)=1.5;
axis(V)

xlabel('N cells X-Y-Z direction');
ylabel('|x_N|_{\inf}/|x_0|')
title('Convergence: 1 F-cycle, constant coef. Laplacian, u=(x^4-x^2) on (0,1)^3, N=32/PE')
grid
print(gcf,'-djpeg100','converg_Fcycles_032N')
%
% plot convergence test - V-256
%
figure; 
load errors_V_256N;
[mv256 n] = size(errors_V_256N);
q2 = ones(mv256+1,1); q1 = q2;

N_V256 = 512*ones(mv256+1,1);
for i=1:mv256
  N_V256(i+1) = 2*N_V256(i);
  s_h = N_V256(i);
  q2(i) = 1/s_h^2*1.e-1;
end

loglog(N_V256(1:mv256),errors_V_256N(:,1),'r*--'),hold on
loglog(N_V256(1:mv256),errors_V_256N(:,2),'go-.'),hold on
loglog(N_V256(1:mv256),q2(1:mv256),'k--'),hold on
%loglog(N_V256(1:mv256),q1(1:mv256),'b-.'),hold on
set(gca,'XTick',N_V256(1:mv256))
legend('Error','Residual','Quadradic',1)
V=axis;
V(1)=N_V256(1);
V(2)=N_V256(mv256);
V(3)=errors_V_256N(mv256,1)/2;
V(4)=errors_V_256N(1,2)*2;
axis(V);

xlabel('N cells X-Y-Z direction');
ylabel('|e_N|_{inf}')
title('Convergence: V-cycles (rtol=1.^{-6}, constant coef. Lapalcian, u=(x^4-x^2) on (0,1)^3, N=256/PE')
grid
print(gcf,'-djpeg100','converg_Vcycles_256N')
%
% finish off error plot
%
%
figure
loglog(N_F256(1:mf256),errors_F_256N(:,1),'bs-.'),hold on
loglog(N_V256(1:mv256),errors_V_256N(:,1),'rx--'),hold on
legend('1 F-cycle w/ V(2,2), N=256/PE, one solve','V-cycles w/ V(2,2), N=256/PE, rtol=1.^{-8}',1),
set(gca,'XTick',N_V256(1:mv256))
V=axis;
V(1)=N_V256(1);
V(2)=N_V256(mv256);
V(3)=errors_V_256N(mv256)*.9;
V(4)=errors_V_256N(1)*1.2;
%axis(V);
xlabel('N cells X-Y-Z direction');
ylabel('|error|_{\inf}')
title('Error convergence: V-cycles rtol=1.^{-6}, constant coef. Laplacian, u=(x^4-x^2) on (0,1)^3, N=256/PE')
grid
print(gcf,'-djpeg100','converg_errors')
%
% times
%
figure
load times_F_256N
load times_F_032N
load times_V_256N
pes = [8 64 512 4096 8*4096 64*4096 512*4096];
[m n] = size(times_F_256N);
semilogx(pes(1:m),times_F_256N,'r+-.'),hold on
semilogx(pes(1:m),times_F_032N,'bs--'),hold on
semilogx(pes(1:m),times_V_256N,'g*-'),hold on
V=axis;
V(1)=pes(1);
V(2)=pes(m);
V(3)=0;
V(4)=times_V_256N(m)*1.2;
axis(V);
set(gca,'XTick',pes(1:m))
legend('1 F-cycle w/ V(2,2), N=256/PE, one solve','1 F-cycle w/ V(2,2), N=32/PE, 512 solves','V-cycles w/ V(2,2), N=256/PE, rtol=1.^{-8}',2),
xlabel('N PEs (Hopper)');
ylabel('Time')
title('Times: V-cycles (rtol=1.^{-6}, constant coef. Laplacian, u=(x^4-x^2) on (0,1)^3, N=256/PE')
grid
print(gcf,'-djpeg100','weak_scaling_hopper')
%
