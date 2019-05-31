function [S_out, N_out] = subject_group_processing(Sig_Pro_type,S,fs_S,N,fs_N,phi,Use_Shadow_filtering,Use_HL_Simulations)




Parser_Cell = strsplit(Sig_Pro_type,'_');


   S_out(:,2)       = process_signal(S(:,2),fs_S,{Parser_Cell{2} Parser_Cell{4}},phi,Use_Shadow_filtering,Use_HL_Simulations,'S');
   temp             = process_signal(S(:,1),fs_S,{Parser_Cell{1} Parser_Cell{3}},phi,Use_Shadow_filtering,Use_HL_Simulations,'S');
   
   S_out(:,1)       = zeropad(temp, length(S_out(:,2)));
   
   N_out(:,2)       = process_signal(N(:,2),fs_N,{Parser_Cell{2} Parser_Cell{4}},phi,Use_Shadow_filtering,Use_HL_Simulations,'N');
   temp             = process_signal(N(:,1),fs_N,{Parser_Cell{1} Parser_Cell{3}},phi,Use_Shadow_filtering,Use_HL_Simulations,'N');
   
   N_out(:,1)       = zeropad(temp, length(N_out(:,2)));

end 

function output = process_signal(input,fs,processing_type,phi,Use_Shadow_filtering,Use_HL_Simulations,Signal_type)
    switch processing_type{1}
    
        case 'NH'
            output  = input;
        case 'CI'
            output  = CI_Sim((input),fs,processing_type{2});
        case 'HL'
            output  = hearing_loss_sim([input, input],fs,processing_type{2},phi,Use_Shadow_filtering,Use_HL_Simulations,Signal_type);
            output  = output(:,1);
        case 'Mon'
            output = zeros(length(input),1);
    end 
end


function output = zeropad(signal,len_new)
    len                 = len_new - length(signal);
    if len_new > length(signal)
        output          = [signal;zeros(len,size(signal,2))];
    else
        output          = signal(1:end+len);
    end 
end
% switch subject_group
%     
% %         if phi==0
% %             S_out     	= audioread('S000_HL1.wav');
% %             N_out      	= audioread('N000_HL1.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL1.wav');
% %             N_out      	= audioread('Np90_HL1.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL1.wav');
% %             N_out      	= audioread('Nm90_HL1.wav');
% %         end
% 
%     case 'NH_bin'
%         S_out       = S; 
%         N_out       = N;
%         
%     case 'NH_mon'
%         S_out(:,1) 	= S(:,1);
%         N_out(:,1)	= N(:,1);
%         
%         S_out(:,2)	= zeros(length(S_out),1);
%         N_out(:,2)	= zeros(length(N_out),1);
%         
%     case 'MED-EL_mon'
%         S_out(:,2)	= CI_Sim(S(1:end-1,2),fs_S,'MED-EL');
%         N_out(:,2)	= CI_Sim(N(1:end-1,2),fs_N,'MED-EL');
%         
%         S_out(:,1)	= zeros(length(S_out),1); % makinging sure both parts of the signal have the same length
%         N_out(:,1)	= zeros(length(N_out),1);
%         
%     case 'MED-EL_SSD'
%         S_out(:,2)	= CI_Sim(S(1:end-1,2),fs_S,'MED-EL');
%         N_out(:,2) 	= CI_Sim(N(1:end-1,2),fs_N,'MED-EL');
%         
%         S_out(:,1)	= zeropad(S(:,1), length(S_out));
%         N_out(:,1)	= zeropad(N(:,1), length(N_out));
%         
%         % Cochlear
%     case 'Cochlear_mon'
%         S_out(:,2)	= CI_Sim(S(1:end-1,2),fs_S,'Cochlear');
%         N_out(:,2)	= CI_Sim(N(1:end-1,2),fs_N,'Cochlear');
%         
%         S_out(:,1)	= zeros(length(S_out),1); % makinging sure both parts of the signal have the same length
%         N_out(:,1)	= zeros(length(N_out),1);
% 
%     case 'Cochlear_SSD'
%         S_out(:,2)	= CI_Sim (S(:,2), fs_S, 'Cochlear');
%         N_out(:,2)	= CI_Sim (N(:,2), fs_N, 'Cochlear');
%         
%         S_out(:,1)	= zeropad(S(:,1), length(S_out));
%         N_out(:,1)	= zeropad(N(:,1), length(N_out));
%         
%         % Hearing Aid 1
%     case 'HG1_mon'
%         
%         S_out     	= hearing_loss_sim (S,fs_S,'HG1');
%         N_out      	= hearing_loss_sim (N,fs_N,'HG1');
% %         if phi==0
% %             S_out     	= audioread('S000_HL1.wav');
% %             N_out      	= audioread('N000_HL1.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL1.wav');
% %             N_out      	= audioread('Np90_HL1.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL1.wav');
% %             N_out      	= audioread('Nm90_HL1.wav');
% %         end
%         S_out(:,2) = zeropad(0,length(S_out(:,1)));
%         N_out(:,2) = zeropad(0,length(N_out(:,1)));
%         
%     case 'HG1_MED-EL'
%         S_out       = hearing_loss_sim(S, fs_S, 'HG1');
%         N_out		= hearing_loss_sim(N, fs_N, 'HG1');
% %         if phi==0
% %             S_out     	= audioread('S000_HL1.wav');
% %             N_out      	= audioread('N000_HL1.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL1.wav');
% %             N_out      	= audioread('Np90_HL1.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL1.wav');
% %             N_out      	= audioread('Nm90_HL1.wav');
% %         end
%         S_out(:,2)	= zeropad(CI_Sim(S(1:end-1,2),fs_S,'MED-EL'),length(S_out));
%         N_out(:,2)	= zeropad(CI_Sim(N(1:end-1,2),fs_N,'MED-EL'),length(N_out));
%         
%     case 'HG1_Cochlear'
%         S_out       = hearing_loss_sim(S, fs_S, 'HG1');
%         N_out		= hearing_loss_sim(N, fs_N, 'HG1');
% %         if phi==0
% %             S_out     	= audioread('S000_HL1.wav');
% %             N_out      	= audioread('N000_HL1.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL1.wav');
% %             N_out      	= audioread('Np90_HL1.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL1.wav');
% %             N_out      	= audioread('Nm90_HL1.wav');
% %         end
%         S_out(:,2)	= zeropad(CI_Sim(S(:,2),fs_S,'Cochlear'),length(S_out));
%         N_out(:,2)	= zeropad(CI_Sim(N(:,2),fs_N,'Cochlear'),length(N_out));
% 
%         % Hearing Aid 2
%     case 'HG2_mon'
%         S_out       = hearing_loss_sim(S, fs_S, 'HG2');
%         N_out		= hearing_loss_sim(N, fs_N, 'HG2');
%             
% %         if phi==0
% %             S_out     	= audioread('S000_HL2.wav');
% %             N_out      	= audioread('N000_HL2.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL2.wav');
% %             N_out      	= audioread('Np90_HL2.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL2.wav');
% %             N_out      	= audioread('Nm90_HL2.wav');
% %         end
%         S_out(:,2)  = zeropad(0,length(S_out(:,1)));
%         N_out(:,2)  = zeropad(0,length(N_out(:,1)));
%         
%     case 'HG2_MED-EL'
%         S_out       = hearing_loss_sim(S, fs_S, 'HG2');
%         N_out		  = hearing_loss_sim(N, fs_N, 'HG2');
% %         if phi==0
% %             S_out     	= audioread('S000_HL2.wav');
% %             N_out      	= audioread('N000_HL2.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL2.wav');
% %             N_out      	= audioread('Np90_HL2.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL2.wav');
% %             N_out      	= audioread('Nm90_HL2.wav');
% %         end
%         S_out(:,2)	= zeropad(CI_Sim(S(1:end-1,2),fs_S,'MED-EL'),length(S_out));
%         N_out(:,2)	= zeropad(CI_Sim(N(1:end-1,2),fs_N,'MED-EL'),length(N_out));
% 
%         
%     case 'HG2_Cochlear'
%         S_out       = hearing_loss_sim(S, fs_S, 'HG2');
%         N_out		= hearing_loss_sim(N, fs_N, 'HG2');
% 
% %         if phi==0
% %             S_out     	= audioread('S000_HL2.wav');
% %             N_out      	= audioread('N000_HL2.wav');
% %         elseif phi==90
% %             S_out     	= audioread('Sp90_HL2.wav');
% %             N_out      	= audioread('Np90_HL2.wav');
% %         else
% %             S_out     	= audioread('Sm90_HL2.wav');
% %             N_out      	= audioread('Nm90_HL2.wav');
% %         end
% 
%         S_out(:,2)	= zeropad(CI_Sim(S(:,2),fs_S,'Cochlear'),length(S_out));
%         N_out(:,2)	= zeropad(CI_Sim(N(:,2),fs_N,'Cochlear'),length(N_out));
% 
%     case 'test'
% %         S_out(:,2)	= CI_Sim(S(1:end-1,2),fs_S,'MED-EL');
% %         N_out(:,2)	= CI_Sim(N(1:end-1,2),fs_N,'MED-EL');
% %         
% %         S_out(:,1)	= CI_Sim(S(1:end-1,1),fs_S,'MED-EL');
% %         N_out(:,1)	= CI_Sim(N(1:end-1,1),fs_N,'MED-EL');
% %         S_out(:,2)    = S(:,2);
% %         N_out(:,2)    = N(:,2);
%         S_out(:,1)      = hearing_loss_sim (S,fs_S,'test');
%         N_out(:,1)      = hearing_loss_sim (N,fs_N,'test');
% 
%         S_out(:,2) = zeropad(0,length(S_out(:,1)));
%         N_out(:,2) = zeropad(0,length(N_out(:,1)));
% 
% end
% end
