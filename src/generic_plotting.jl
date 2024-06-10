# Establish our plotting functions with dummy methods that will never be used
# (using zero-arg methods for this isn't sufficient since concordialine 
#  and concordialine! may have real zero-arg methods)
struct NotUsed end

function concordiacurve(x::NotUsed) end

function concordiacurve!(x::NotUsed) end

function concordialine(x::NotUsed) end

function concordialine!(x::NotUsed) end

function rankorder(x::NotUsed) end

function rankorder!(x::NotUsed) end