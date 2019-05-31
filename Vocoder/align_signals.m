function [aligned_signal,xlags_in_out] = align_signals(shifted_signal, ref_signal)
        shifted_signal = shifted_signal(:); %Make vector
        ref_signal = ref_signal(:); %Make vector
        [xcorrelation, xlags] = xcorr(real(shifted_signal),real(ref_signal));
        [~, xInd] = max(xcorrelation);
        xlags_in_out = xlags(xInd);
        % adjust ref signal accordingly
            if xlags_in_out < 0
                aligned_signal = [ref_signal(abs(xlags_in_out)+1:end); zeros(abs(xlags_in_out),1)];
            elseif xlags_in_out > 0
                aligned_signal = [zeros(abs(xlags_in_out),1); ref_signal(1:end-abs(xlags_in_out))];
            elseif xlags_in_out == 0
                aligned_signal = ref_signal;
            end
        aligned_signal = aligned_signal(:); %make vector
end