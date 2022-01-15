function [fidexmult,X_IT_new,X_IT_old] = Multi_Sensor_FASTEST_IF_Recover_ICCD(Sig,N_S,win_length, num, delta,L,thr,Thr,step,FFT_len,iiii)
% SET OF FRACTIONAL WINDOWS
w=gausswin(win_length,1);
l=0:win_length-1;
Siggg=Sig;
X_IT=Sig;
for iii=1:FFT_len
    WW(iii,:)=exp(-1i*(iii)*2*pi*l/FFT_len);
end
e_max=0;
%tic;
i=0;
window_rot=zeros(2*L+1,win_length);
for k=-L+1:1:L-1
    i=i+1;
    window_rot(i,:)=frft(w,0.85* k/L);%0.05
end
for jj=1:N_S
    indexs{jj}=(find(Siggg(jj,:)==0));
end
%save('window_rot','window_rot');
%load('window_rot');
%toc;
for tttt=1:1
    for jj=1:N_S
        X_IT(jj,indexs{jj})=0;
        
    end
    for iii=1:num
        for jj=1:N_S
            Sig_extended(jj,:)=[zeros(1,floor(win_length/2)) Sig(jj,:) zeros(1,floor(win_length/2))];
        end
        if N_S>1
        Siga=filter(ones(1,win_length),1,sum(abs(Sig)));
        else
                    Siga=filter(ones(1,win_length),1,(abs(Sig)));

        end
        % Siga=
        
        [~,t_start]=max(Siga(floor(win_length/2)+1:end-floor(win_length/2)));
        
        t_start=t_start(1)+floor(win_length/2);
        % t_start=46;
        v_m=0;
        for i=1:2*L+1
            FF=zeros(1,length(Sig));
            for jj=1:N_S
                FF=FF+abs(fft(Sig(jj,t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),length(Sig)));
            end
            [v,index]=max(FF(1:end/2));
            if v_m<v
                freq_start=index-1;  %peak_frequency location
                frac_wind_index_start=i;   %chirp rate location
                v_m=v;
            end
        end
        v_oldd=v_m;
        
        %t_start=t_start-floor(win_length/2)-1;
        IF=zeros(1,length(Sig))-1;
        IF(t_start)=freq_start;
        
        
        for it=1:2
            f0=freq_start;
            frac_wind_index=frac_wind_index_start;
            t_index=t_start;
            while and(t_index>1,t_index<length(Sig))
                if it==1
                    
                    if t_index+step>length(Sig)
                        %f0=f0+(length(Sig)/win_length)*(length(Sig)-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                        f0=f0+(FFT_len/win_length)*(length(Sig)-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                        t_index=length(Sig);
                        %             IF(t_index)=round(f0);
                        
                        %  break;
                    else
                        t_index=t_index+step;
                        f0=f0+(FFT_len/win_length)*step*tan(1*pi*(frac_wind_index-(L))/(2*L));
                    end
                    
                else
                    if t_index-step<1
                        f0=f0-(FFT_len/win_length)*abs(1-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                        t_index=1;
                        %            IF(t_index)=round(f0);
                        %break;
                        
                    else
                        t_index=t_index-step;
                        f0=f0-(FFT_len/win_length)*step*tan(1*pi*(frac_wind_index-(L))/(2*L));
                    end
                    
                    
                end
                f0=round(f0);
                
                v_m=0;
                k=f0-delta:1:f0+delta;
                k(k>FFT_len/2)=FFT_len/2;
                k(k<=0)=1;
                if frac_wind_index<2
                    frac_wind_index=2;
                elseif frac_wind_index>2*L-1
                    
                    frac_wind_index=2*L-1;
                end
                for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS
                    
                    for jj=1:N_S
                        w_signal(jj,1:win_length)=(Sig_extended(jj,t_index:t_index+win_length-1).*window_rot(i,:)); %WINDOWED SIGNAL
                    end
                    for jjj=1:length(k)%jj
                        VVV=0;
                        for jj=1:N_S
                            VVV=VVV+(abs(sum(w_signal(jj,:).*WW(k(jjj),:))));
                        end
                        V(jjj)=VVV;
                    end
                    [v,index]=max(V);
                    if v_m<v
                        f0=k(index);  %peak_frequency location
                        frac_wind_index=i;   %chirp rate location
                        v_m=v;
                    end
                end
                if v_m<Thr*v_oldd
                    
                    break;
                end
                
                IF(t_index)=f0;
            end
        end
        %IF(t_index)=f0;
        
        ind=find(IF>0);
        IF=interp1(ind,IF(ind),1:length(Sig));
        IF(isnan(IF))=0;
        IF=IF/(1*FFT_len);
        
        
        
        %figure; plot(IF)
        Phase=2*pi*filter(1,[1 -1],IF);
        s_dechirp=exp(-1i*Phase);
        
        
        LL=delta;
        %TF filtering for each sensor
        fidexmult(iii,:) = IF;
        ssss=0;
        
        for jj=1:N_S
            
          [s11,~] = ICCD_sparse(Sig(jj,:),1,IF,1,3,0,iiii{jj});

            
           
            if tttt==1
                Sig(jj,iiii{jj})=Sig(jj,iiii{jj})-s11(iiii{jj});
            else
                Sig(jj,:)=Sig(jj,:)-s11;
                
            end
            
            XAA(jj,iii,:)=s11;
            X_IT(jj,indexs{jj})=X_IT(jj,indexs{jj})+s11(indexs{jj});
            % end
        end
       % I=HTFD_new1(s11,3,8,64);
        %figure;imagesc(I)
        if ssss<e_max*thr
            break;
        end
        if ssss>e_max
            e_max=ssss;
        end
        
        
    end
    for jj=1:N_S
    [XAB(jj,:,:),~,~] = ICCD_sparse(Siggg(jj,:),1,fidexmult,1,3,0,iiii{jj});
    end

     XAA=XAB;
 for i=1:1
    for kkk1=1:N_S
        for jjj1=1:num
            B=reshape(XAA(kkk1,jjj1,:),1,128);
            A=reshape(XAA(1,jjj1,:),1,128);
              H(kkk1,jjj1)= A*B'/(B*B');
        end
    end
    %abs(H)
        [XAA,~] = ICCD_sparse_multi(Siggg,1,fidexmult,3,0,iiii,conj(H));

 end
        
for j=1:N_S
     X_IT_new(j,:)=sum(XAA(j,:,:));
   %  X_IT_new(j,iiii{j})=Siggg(j,iiii{j});
     X_IT_old(j,:)=sum(XAB(j,:,:));
    %      X_IT_old(j,iiii{j})=Siggg(j,iiii{j});

end

   Sig=X_IT_new; 
    %Sig=Sig_out1;
end
