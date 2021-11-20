%====================================================================
%                          MIAA Project
%====================================================================
%                                          Original Coder: Y. Chen
%                                             Modified by: J. Ling
%                                                    Date: 08/12/2008
%--------------------------------------------------------------------
%   First thing first. This simple version is to get familiar with 
%   the methods proposed in the paper.
%--------------------------------------------------------------------

    
   
    function [Sig]=recover_component(dataGiven,GivenIX,MissingIX)
    % number of the missing samples (out of 100)
%    Num_Missing_Samples = 5;
    % number of the Monte-Carlo trials
    % frequency locations of the true line spectra
    
    Fs       = 1;           % sampling frequency
    T        = 1/Fs;        % sampling period
    Fre_Grid = 1000;        % frequency grid $K$ in the paper, so the finest grid is 1/1000 Hz 
    
    
%        iMonteCarlo
        %signal generation
  %      noise = wgn(1,Data_Len,-20,'complex');
 %       truePhase = unifrnd(-pi,pi,1,4);
    %   truePhase = zeros(1,4);
%        data = trueAmp.*exp(j*truePhase)*exp(2j*pi*trueFreq'*Tem_Grid) + noise;
            
            
            
            %%--------------------------------------------------------
            %%The IAA spectral estimation from from the available data
            %%--------------------------------------------------------
            %figure; plot(abs(Amplitude_Spectrum_IAA));
            Fre_Line  = 2*pi*(1:Fre_Grid)/Fre_Grid;
            A_m  =  exp(1i*MissingIX'*T*Fre_Line);
        
            % apply MIAA-1 and MIAA-2 methods to recover the missing data
            % see Yubo's letter for details
           %Sig = miaa1(Amplitude_Spectrum_IAA, Time_Line, Fre_Line);        
            %Sig = miaa2(dataGiven, A_g, A_m, Amplitude_Spectrum_IAA, Fre_Line);
             %Sig = miaa3(dataGiven, GivenIX, Fre_Grid, 15, A_m);
                          Sig = miaa3(dataGiven, GivenIX, Fre_Grid, 5, A_m);

          %Sig=ym_miaa1;
            
          %  MSE1(IXIX) = MSE1(IXIX) + norm(dataMissing - ym_miaa1)^2/Num_Missing_Samples;   % MSE for MIAA-1
          %  MSE2(IXIX) = MSE2(IXIX) + norm(dataMissing - ym_miaa2)^2/Num_Missing_Samples;   % MSE for MIAA-2
          %  MSE3(IXIX) = MSE3(IXIX) + norm(dataMissing - ym_miaa3)^2/Num_Missing_Samples;   % MSE for MIAA-3
          %  MSE4(IXIX) = MSE4(IXIX) + norm(dataMissing - ym_mapes)^2/Num_Missing_Samples;   % MSE for MAPES
        
%        title = ['LOG\PACKAGE_' num2str(iMonteCarlo)];
%        save(title, 'IAA_T', 'TIME1', 'TIME2', 'TIME3', 'TIME4'); %, 'MSE1', 'MSE2', 'MSE3', 'MSE4');
    end
    
    
