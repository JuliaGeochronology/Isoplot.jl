module MakieExt

    using Isoplot
    using Makie
    using Measurements
    import Isoplot: concordialine, concordialine!, concordiacurve, concordiacurve!

    Makie.convert_arguments(P::Type{<:Poly},e::Ellipse) = convert_arguments(P,as_points(e.x,e.y))
    Makie.convert_arguments(P::Type{<:Poly},a::Analysis) = convert_arguments(P,Ellipse(a))
    Makie.plottype(::Ellipse) = Poly
    Makie.plottype(::Analysis) = Poly

    



    function as_points(v1,v2)
        
        if length(v1) == length(v2)
            pointList = Point2f[]
            for i in 1:lastindex(v1)
                push!(pointList,Point2f(v1[i],v2[i]))
            end
            return pointList
        else
            Throw(DimensionMismatch("Arguments must be the same length"))
        end
    end

    # Plot a line between two times in Wetherill Concordia space
    @recipe(ConcordiaLine, t₀, t₁) do scene
        Attributes(
            concordiaType = :Wetherril,
            color = :black,
            linewidth = 1,
            errorcolor = (:black,0.15)
            
            

        )
    end

    function Makie.plot!(cline::ConcordiaLine)
        r75₀ = ratio(cline.t₀[], λ235U)
        r68₀ = ratio(cline.t₀[], λ238U)
        r75₁ = ratio(cline.t₁[], λ235U)
        r68₁ = ratio(cline.t₁[], λ238U)
    
        slope = 0.0 ± 0.0
        intercept = 0.0 ± 0.0
        if cline[:concordiaType][] == :Wetherril
            slope = (r68₁-r68₀)/(r75₁-r75₀)
            
            intercept = r68₀ - r75₀*slope
        elseif cline[:concordiaType][] == :TeraWasserburg
            #TODO
        else   
            throw(ArgumentError("concordiaType must be :Wetherril or :TeraWasserburg"))
        end
        xmin = val(r75₀)
        xmax =val(r75₁)
        x = Observable(collect(range(xmin, xmax, length=50)))
        y = intercept .+ slope .* x[]
        y_val = Observable(val.(y))
        y_err = err.(y)
        
        upperError = as_points(x[],y_val[] .+ y_err)
        lowerError = as_points(reverse(x[]),reverse(y_val[] .- y_err))
        errorPoly = Observable([upperError;lowerError])
        poly!(cline,errorPoly,color=cline[:errorcolor][],strokewidth =0)
        lines!(cline,x,y_val,color = cline[:color][], linewidth = cline[:linewidth][])
        return cline
    end

    @recipe(ConcordiaCurve,t₀, t₁) do scene
        Attributes(
            concordiaType = :Wetherril,
            color = :black,
            linewidth = 1,
            errorcolor = (:black,0.15),
            agefontsize = 10,
            ageticksize = 5,
            agetickcolor = :black
        )
    end


    function Makie.plot!(ccurve::ConcordiaCurve)
        # Uncertainty of 235 decay constant relative to the 238 decay constant
        σₜ = λ235U_jaffey.val .* sqrt((λ238U.err/λ238U.val).^2 + (λ235U_jaffey.err/λ235U_jaffey.val).^2) # 1/Years

        # Plot the concordia curve
        tlim = [ccurve.t₀[],ccurve.t₁[]]
        
        # xl = [ccurve.t₀[],ccurve.t₁[]]
        # xl, yl = Plots.xlims(hdl), Plots.ylims(hdl) # Note current size of figure
        # tlim = age.(max.(xl, 0.0), λ235U_jaffey.val) # Calculate time range of current window
        
        dt = tlim[2] - tlim[1]
        tmin = max(tlim[1]-0.1dt, 0.0)
        tmax = tlim[2]+0.1dt
        t = range(tmin, tmax, length=1000) # Time vector, including padding

        xratio = Observable(Float64[])
        yratio = Observable(Float64[])
        errx = Float64[]
        erry = Float64[]

        scale = floor(log10(tlim[2]-tlim[1])) # Characteristic timescale (order of magnitude)
        trange = round.(tlim./10.0^scale) # Minimum and maximum time to a round number
        majorstep = 0.5
        tticks = (trange[1]:majorstep:trange[2]).*10.0^scale # Time ticks, to a round number
        tickx = Observable(Float64[])
        ticky = Observable(Float64[])


        if ccurve[:concordiaType][] == :Wetherril
            xratio[]= ratio.(t, λ235U_jaffey.val) # X axis values
            yratio[] = ratio.(t, λ238U.val)# Y axis values

            errx = [ratio.(t, λ235U_jaffey.val-σₜ*2); reverse(ratio.(t, λ235U_jaffey.val+σₜ*2))]
            erry = [yratio[]; reverse(yratio[])]

            tickx[] = ratio.(tticks, λ235U_jaffey.val) # X axis values
            ticky[] = ratio.(tticks, λ238U.val)# Y axis values
            
        elseif ccurve[:concordiaType][] == :TeraWasserburg
            xratio[]= ratio.(t, λ235U_jaffey.val) # X axis values
            yratio[] = ratio.(t, λ238U.val)# Y axis values

            errx = [ratio.(t, λ235U_jaffey.val-σₜ*2); reverse(ratio.(t, λ235U_jaffey.val+σₜ*2))]
            erry = [yratio[]; reverse(yratio[])]

            tickx[] = ratio.(tticks, λ235U_jaffey.val) # X axis values
            ticky[] = ratio.(tticks, λ238U.val)# Y axis values
        else   
            throw(ArgumentError("concordiaType must be :Wetherril or :TeraWasserburg"))
        end
        errorPts = Observable(as_points(errx,erry))
        poly!(ccurve,errorPts,color=ccurve[:errorcolor][],strokewidth =0)
        lines!(ccurve,xratio, yratio,color = ccurve[:color][], linewidth = ccurve[:linewidth][])
        
        # # Plot age markers with text labels
        scatter!(ccurve,tickx,ticky,markersize = ccurve[:ageticksize][],color = ccurve[:agetickcolor][],transform_marker=true)
        xoffset = (maximum(tickx[])-minimum(tickx[]))/200
        yoffset = (maximum(ticky[])-minimum(ticky[]))/100
    
        t = (minimum(xratio[])+8*xoffset) .< tickx[] .< maximum(xratio[])
        ticklabels = Observable(string.(round.(Int, tticks[t])))
        #can probably make these mobile based on zoom level
        tickLabelX = Observable(tickx[][t].-xoffset)
        tickLabelY = Observable(ticky[][t].+yoffset)
        text!(ccurve,tickLabelX,tickLabelY,text=ticklabels,fontsize = ccurve[:agefontsize][],transform_marker=true)
        # Plots.annotate!(hdl, r75tticks[t].-xoffset,r68tticks[t].+yoffset,ticklabels)

        return ccurve
    end
end


