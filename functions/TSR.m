function sequence = TSR(TS_coeff,t,n,method)

% function [FT,FD,SD] = TSR(TS_coeff,t,n,method)

    N = size(TS_coeff,1);
    l = size(TS_coeff,2);
    n_f = length(t);
    sequence  = zeros(N,l,n_f,'single'); 
%     FT  = zeros(N,l,n_f); 
%     FD  = zeros(N,l,n_f); 
%     SD  = zeros(N,l,n_f); 
    
    switch method
        case 'FT' % Fitted thermal sequence
            for f = 1:n_f % Frames Loop
               log_tn = log(t(f));
               for k = 1:n+1
                     sequence(:,:,f) = sequence(:,:,f) + TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-1);
               end
            end
            sequence = exp(sequence);
        case 'FD' % first derivative of FT
            for f = 1:n_f % Frames Loop
               log_tn = log(t(f));
               for k = 1:n+1
                    if k>1
                        sequence(:,:,f) = sequence(:,:,f) + (k-1)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-2);
                    end
               end
            end
            sequence = exp(sequence);
        case 'SD' % second derivative of FT
            for f = 1:n_f % Frames Loop
               log_tn = log(t(f));
               for k = 1:n+1
                    if k>2 
                        sequence(:,:,f) = sequence(:,:,f) + (k-1)*(k-2)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-3);
                    end
               end
            end
            sequence = exp(sequence);
        case 'ALL' %all cases
            sequenceD  = zeros(N,l,n_f,'single');
            sequenceDD  = zeros(N,l,n_f,'single');
            for f = 1:n_f % Frames Loop
                log_tn = log(t(f));
                for k = 1:n+1
                    sequence(:,:,f) = (sequence(:,:,f) + TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-1)); %FT
                    if k >1 sequenceD(:,:,f) = (sequenceD(:,:,f) + (k-1)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-2)); end %FD
                    if k >2 sequenceDD(:,:,f) = (sequenceDD(:,:,f) + (k-1)*(k-2)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-3)); end%SD
                end
            end
            sequence = {exp(sequence), sequenceD, sequenceDD}; %return as a cell
    end
%     for f = 1:n_f % Frames Loop
%         log_tn = log(t(f));
% 
%         for k = 1:n+1
%             
%              if k>1
%                 FD(:,:,f) = FD(:,:,f) + (k-1)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-2);
%              end
%              if k>2 
%                 SD(:,:,f) = SD(:,:,f) + (k-1)*(k-2)*TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-3);
%              end
%              FT(:,:,f) =    FT(:,:,f) + TS_coeff(:,:,-k+(n+2))*(log_tn)^(k-1);
%         end
% 
%     end



end



