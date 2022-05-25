function [Precision, Recall, JS, SA, Score, k]=pre_rec(A,B)
[l1, l2] = size(A);
A=A(:);
B=B(:);
C=A(:)+B(:);
TN=sum(C==0);
TP=sum(B(A==1)==1); % AB %true positive
FN=sum(B==1)-TP; % true negative
FP=sum(A==1)-TP; % false positvie
TT=l1*l2-sum(C==0);
P0=(TP+TN)/(TP+FN+FP+TN);
a1=TP+FN;
a2=FP+TN;
b1=TP+FP;
b2=TN+FN;
Pe=(a1*b1+a2*b2)/((TP+FN+FP+TN)*(TP+FN+FP+TN));
k=(P0-Pe)/(1-Pe);
SA=P0;
Precision=TP/(TP+FP);
Recall=TP/(TP+FN);
JS=TP/TT;
Score=2*Precision*Recall/(Recall+Precision);
return