function [ry cx mx]=cmbn9(mx0)
global fs fgn lgn scl cmap sz rot ump pc brn

[mx r c]=mx2d(mx0, sz);

ry=ump.*(r-1)./1000; cx=ump.*(c-1)./1000;

imf(brn); colormap gray; title('mx-amp');

nm=strvcat('S0','S1','S2','S3','S4','S5','S6','S7','S8')

for k0=1:length(r)
    hold on; text(cx(k0), ry(k0), cutstr(nm(k0,:)), 'Color' , 'red','HorizontalAlignment', 'center')
end




% 
% [ds_mx rc_mx]=rtmx(mx0); title('Amplitude')%% retinotopic map of max amplitude va
% [ds_rx rc_rx]=rtmx(rx0); title('Ratio') %% retinotopic map of ratio amplitude values 
% L_mx=kntr(mx0, kn, rc_mx); title('Amplitude') %% 
% L_rx=kntr(rx0, kn, rc_rx); title('Ratio') %%
% 
% imf(brn); colormap gray; title('mx-amp');
% 
% cl=[1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 0; 1 1 1; 0 0 0]; %; 0 0 0]; %0.6350, 0.0780, 0.1840];
% ln{1}=[4 6 7]; ln{2}=[4 5 6 7 8]; ln{3}=[5 6 8]; ln{4}=[6 7]; ln{5}=8; ln{6}=8; 
% for k0=1:size(ds_mx,1);
%     %     hold on; plot(rc_mx(k0,2), rc_mx(k0,1), '*k'); hold on; plot(rc_mx(k0,2), rc_mx(k0,1), 'dw')
%     hold on; plot(rc_mx(k0,2), rc_mx(k0,1), 'Marker', '*', 'Color', cl(k0,:));
%     hold on; plot(rc_mx(k0,2), rc_mx(k0,1), 'Marker', 'd', 'Color', cl(k0,:));
%     %if k0<=length(ln)   %%exist(['ln{' num2str(k0) '}'])
%     if k0<=length(ln)
%         for k1=1:length(ln{k0})
%             hold on; plot([rc_mx(k0,2) rc_mx(ln{k0}(k1),2)], [rc_mx(k0,1) rc_mx(ln{k0}(k1),1)], 'k');
%         end
%     end
% end
% 
% 
% ds{1}=ds_mx; ds{2}=ds_rx;
% rf{1}=rc_mx(:,1); rf{2}=rc_rx(1,:);
% cf{1}=rc_mx(:,2); cf{2}=rc_rx(:,2);
% L{1}=L_mx; L{2}=L_rx;
