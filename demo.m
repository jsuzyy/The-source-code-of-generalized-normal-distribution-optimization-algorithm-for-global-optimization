clc
clear all

lb = -100.*ones(1,10);
ub = 100.*ones(1,10);
maxit = 100;
objf= @Sphere;

n = 30;
d = 10;

[BestCost,BestValue,Best]=GNDO(objf,n,d,lb,ub,maxit);
plot(BestCost,'r','linewidth',2)
xlabel('The number of iterations','Fontname','Times New Roma','fontsize',15);
ylabel('Fitness value','Fontname','Times New Roman','fontsize',15);