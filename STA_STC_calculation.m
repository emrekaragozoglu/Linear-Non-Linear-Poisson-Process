% history=32;
% thresh=0;
% for i=1:1:18
%     spikes = find(est_data(i).resp>thresh);
%     spikes=spikes.*32;
%     for n=1:1:length(spikes)
%         if spikes(n)>history-1
%             spike_trig_array(n,:)=(est_data(i).resp(spikes(n)./32)/0.1)*transpose(est_data(i).stim(spikes(n)-(history-1):spikes(n)));
%         else
%             num=history-spikes(n);
%             z=zeros(1,num);
%             spike_trig_array(n,:)=(est_data(i).resp(spikes(n)./32)/0.1)*horzcat(z,transpose(est_data(i).stim(1:spikes(n))));
%         end
%     end
%     [row,col] = size(spike_trig_array);
%     STA=sum(spike_trig_array,1)./row;
%     STAmatrix(i,:)=STA;
%     STA=zeros(1,history);
% end
% finalSTA=sum(STAmatrix,1)/.18;
% figure
% plot(finalSTA)


% history=32;
% thresh=0;
% i=1;
% 
% for s=1:1:10
%     spikes=find(transpose(est_data(1).resp_raw(:,s))>thresh);
%     spikes=spikes.*32;
%     for n=1:1:length(spikes)
%         if spikes(n)>history-1
%             spike_trig_array(n,:)=transpose(est_data(1).stim(spikes(n)-(history-1):spikes(n)));
%         else
%             num=history-spikes(n);
%             z=zeros(1,num);
%             spike_trig_array(n,:)=horzcat(z,transpose(est_data(1).stim(1:spikes(n))));
%         end
%     end
%     [row,col] = size(spike_trig_array);
%     STA(s,:)=sum(spike_trig_array,1)./row;
%     STA(s,:)=STA(s,:)./norm(STA(s,:));
% end
% 
% %Calculating the STC
% 
% for s=1:1:10
%     proj=conv(fliplr(STA(s,:)),est_data(1).stim);
%     projected=proj(1:length(proj)-length(STA(s,:)));
%     
%     for i=1:1:length(projected)
%         if projected(i)>0
%             projected(i)=projected(i);
%         else
%             projected(i)=0;
%         end
%     end
%     poiss(s,:)=poissrnd(projected);
% end
% finalpoiss=sum(poiss,1)./10;
% 
% c=0;
% for i=1:32:length(finalpoiss)
%     c=c+1;
%     train(c)=finalpoiss(i);
% end
% 
% figure
% subplot(2,1,1)
% plot(est_data(1).resp)
% subplot(2,1,2)
% plot(train,'r');
% length(find(est_data(1).resp~=0))
% length(find(train~=0))

%Concatination
respconcat=[];
stimconcat=[];
for t=1:1:18
    respconcat=[respconcat transpose(est_data(t).resp)];
    stimconcat=[stimconcat transpose(est_data(t).stim)];
end


history=32;
thresh=0;
i=1;

%Calculating STA
spikes=find(respconcat>thresh);
spikes=spikes.*32;
for n=1:1:length(spikes)
    if spikes(n)>history-1
        spike_trig_array(n,:)=(stimconcat(spikes(n)-(history-1):spikes(n)));
    else
        num=history-spikes(n);
        z=zeros(1,num);
        spike_trig_array(n,:)=horzcat(z,(stimconcat(1:spikes(n))));
    end
end
[row,col] = size(spike_trig_array);
STA=sum(spike_trig_array,1)./row;
STA=STA./norm(STA);
plot(STA)
%%
%Calculating the STC
num_spike = zeros(1,length(spikes));
for i = 1:length(spikes)
    if(spikes(i)/32 > 31)
        num_spike(i) = sum(est_data(1).resp((spikes(i)/32-31):(spikes(i)/32)))*10;  
    end
end

STC = zeros(32,32);
for j=1:1:length(spikes)
    a = spikes(j);
    if a > 31
        STC=STC+(num_spike(j))*transpose((transpose(est_data(1).stim(a-31:a))- STA))*(transpose(est_data(1).stim(a-31:a))-STA);
    end
end
STC=STC./(sum(num_spike));
%%
proj=conv(fliplr(STA),stimconcat);
projected=proj(1:length(proj)-length(STA)+1);

% for i=1:1:length(projected)
%     if projected(i)>0
%         projected(i)=projected(i);
%     else
%         projected(i)=0;
%     end
% end
% poiss=poissrnd(projected);
% 
% finalpoiss=sum(poiss,1)./10;
% 
% c=0;
% for i=1:32:length(projected)
%     c=c+1;
%     train(c)=projected(i);
% end
% 
% figure
% subplot(2,1,1)
% plot(est_data(1).resp)
% subplot(2,1,2)
% plot(train,'r');
% length(find(est_data(1).resp~=0))
% length(find(train~=0))
% figure
c=0;
for i=1:32:length(projected)
    c=c+1;
    decimated(c)=projected(i);
end
figure
scatter(respconcat(1:35040),decimated);

%%
x = [0 0.1 0.2 0.3 0.4];
y = [0.1989 0.1395 0.08826 0.03697 0.01022];



%%

input=randn(1,10000);
mu1=18; 
mu2=12; 
sigma1=1;
sigma2=5;

x=1:20;
gauss1=exp(-(x-mu1).^2/2/sigma1^2);
gauss2=0.2*exp(-(x-mu2).^2/2/sigma2^2);
k=gauss1-gauss2;
k=k./norm(k); %normalize K

%convolve linear filter with input
linOut=conv(input,fliplr(k));
linOut(length(linOut)-18:length(linOut))=[];

instRate=0.1*log(1 + exp(10*linOut));

train=poissrnd(instRate);
train(isnan(train)) = 0;

plot(instRate)
%%

% figure
% subplot(3,1,1)
% plot(est_data(1).stim)
% subplot(3,1,2)
% plot(projected,'r')
% subplot(3,1,3)
% plot(est_data(1).resp)

projectedsq=0.7*log(1+exp(0.6*projected));
%projectedsq=projected.^2;
train=poissrnd(projectedsq);

figure
subplot(2,1,1)
plot(est_data(3).resp)
subplot(2,1,2)
plot(train)

% figure
% train(isnan(train))=0;
% plot(train)

% for i=1:1:length(train)
%     if train(i)>=2
%         traintresh(i)=1;
%     else
%         traintresh(i)=0;
%     end
% end
% 
% figure
% subplot(2,1,1)
% plot(est_data(1).resp)
% subplot(2,1,2)
% plot(traintresh)
