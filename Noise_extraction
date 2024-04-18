function Y_N=Noise_extraction(p11,aa1,YY1,N11,p_r11)
k_pp=1;
if -p11+1+p_r11*k_pp>=p11+1
k_pp=k_pp;
else k_pp=k_pp+1;
end
peri=p11+1;
for i5=1:p_r11*N11-p11
if peri==p_r11+1 peri=1;
end
for j3=1:p11+1
Q_Mat(i5,i5+j3-1)=-aa1(p11+1-j3+1,peri);
end
peri=peri+1;
end
Q_Mat1(:,:)=Q_Mat(:,p11+1:end);
NoiSe1(p11+1:p_r11*N11,1)=Q_Mat1*YY1(p11+1:p_r11*N11,1);
for i2=1:p11
aa_f1=flip(aa1(2:p11+1,i2));
NoiSe1(i2,1)=-dot(YY1(k_pp*p_r11-p11+i2:k_pp*p_r11-1+i2),aa_f1)+YY1(i2,1);
end
Y_N=NoiSe1;
end
