# Here I compute the quantum anomalous hall (B) independent component.
# Take into consideration that you'll need to break C2x in order to do so.
# C2x is σxτx so a sigmaz will break it given rise to a QAH of course, you
# need chern bands to start with

function qah_contribution(a, b)
    _warncrossed(a,b)
    return 


end

function qah_contribution_k(a,b)


end

function _warncrossed(a,b)
    if a == b
        return @code_error
    elseif a == :z || b == :z
        return @code_error
    else nothing end
end
