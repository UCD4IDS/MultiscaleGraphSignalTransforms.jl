using Plots
"""
    GraphSig_Plot(G::GraphSig)

Display a plot of the data in a GraphSig object

### Input Argument
* `G::GraphSig`: an input GraphSig object

"""
function GraphSig_Plot(G::GraphSig)

    # only run this for 1D signals
    fcols = Base.size(G.f, 2)
    if fcols == 0
        G.f = zeros(G.length, 1)
    elseif fcols != 1
        warn("Cannot handle multiple signals at this time!")
        return
    end

    ## Case 1: the graph does not have spatial coordinates or it is too large
    if G.length > 10^4 || G.dim > 3 || G.dim < 1
        if G.length < 10^6
            spy(G.W, show = true)            # not too clear if we can change the dot size.
            scatter!((1:G.length)', (1:G.length)', marker=(:circle, 4), zcolor = G.f, title = G.name)
        else
            warn("Cannot handle a graph with more than 10^6 nodes at this time!")
        end
        return
    end


    ## Case 2: the graph does have spatial coordinates (G.dim = 1, 2, or 3).

    # extract the plot specifications
    symmetric,markercolor,markerstrokealpha,markervaluevaries,notitle,nocolorbar,stemplot,cLim,cmin,cmax,ptsize,linewide,linecolor,verbatim,verbtext,sortnodes = ExtractPlotSpecs(G);
    # set(0, 'DefaultFigureVisible', 'on');

    #################################
    ############# 1-D ###############
    #################################
    if G.dim == 1
        if G.length < 65 || stemplot
            # stem(G.xy, G.f, '-o','LineWidth',2,'Color',linecolor,'MarkerFaceColor',linecolor,'MarkerSize',ptsize(1));
            stem(G.xy, G.f, linetype = [:sticks :scatter], linewidth = 2, linecolor = linecolor, markercolor = linecolor, markersize = ptsize[0])
            # plot(1:G.length , G.f, line = :sticks)
            # plot!(1:G.length , G.f, line = :scatter)
        else
            # plot(G.xy, G.f, '-','Color',linecolor,'LineWidth',linewide);
            plot(G.xy, G.f, linetype = :line, linecolor = linecolor, linewidth = linewide);
        end

        xmin = minimum(G.xy);
        xmax = maximum(G.xy);

        #     dx = xmax-xmin;
        #     if dx < 10*eps
        #         dx = 1;
        #     end
        #     xmin = xmin - 0.1*dx;
        #     xmax = xmax + 0.1*dx;

        ymin = minimum(G.f);
        ymax = maximum(G.f);
        dy = ymax-ymin;
        if dy < 10*eps
            dy = 1;
        end
        ymin = ymin - 0.1*dy;
        ymax = ymax + 0.1*dy;

        # symmetric?
        if symmetric
            ymax = 1.1*max(ymax, -ymin);
            ymin = -ymax;
        end

        #axis([xmin, xmax, ymin, ymax]);
        xaxis!(xlim = (xmin, xmax))
        yaxis!(ylim = (ymin, ymax))

        if cLim
            #axis([xmin, xmax, cmin, cmax]);
            xaxis!(xlim = (xmin, xmax))
            yaxis!(ylim = (cmin, cmax))
        else
            #axis([xmin, xmax, ymin, ymax]);
            xaxis!(xlim = (xmin, xmax))
            yaxis!(ylim = (ymin, ymax))
        end

        # title?
        if ~isempty(G.name) && ~notitle
            title!(G.name);
        end

        # verbatim? (need to update later)
        #if verbatim
        #    eval(verbtext);
        #end


        #################################
        ######### 2-D  or 3-D ###########
        #################################

    else

        #fig = figure('visible','on');

        if sortnodes
            IX = sortperm(abs.(G.f[:,1]),rev = true);
            G.W = G.W[IX,IX];
            G.xy = G.xy[IX,:];
            G.f = G.f[IX,:];
        end

        # plot the graph
        if isempty(linewide)
            gplot(G.W, G.xy, style = :solid, color = linecolor);
        else
            gplot(G.W, G.xy, style = :solid, color = linecolor, width = linewide);
        end
        #hold on

        # plot the nodes

        if ptsize[1] != 0
            # if we're using various point sizes, determine the vector of sizes
            if length(ptsize) > 1#~isscalar(ptsize)
                if cLim
                    ptsize = [ptsize[1] .+ ptsize[2].*abs(G.f)./maximum(abs([cmin,cmax]))];
                else
                    ptsize = [ptsize[1] .+ ptsize[2].*abs(G.f)./maximum(abs(G.f))];
                end
            end

            # plot the nodes for a 2-D graph
            if G.dim == 2
                if ~markervaluevaries
                    scatter!(G.xy[:,1], G.xy[:,2], markersize = ptsize[1], markercolor = markercolor,markerstrokealpha = markerstrokealpha);
                else
                    scatter!(G.xy[:,1], G.xy[:,2], markersize = ptsize[1], zcolor = G.f, markercolor = markercolor,markerstrokealpha = markerstrokealpha);
                end

                # plot the nodes for a 3-D graph
            else
                if ~markervaluevaries
                    scatter!(G.xy[:,1], G.xy[:,2], G.xy[:,3], markersize = ptsize[1], markercolor = markercolor,markerstrokealpha = markerstrokealpha);
                else
                    scatter!(G.xy[:,1], G.xy[:,2], G.xy[:,3], markersize = ptsize[1], zcolor = G.f, markercolor = markercolor,markerstrokealpha = markerstrokealpha);
                end
            end

            # colorbar?
            #if !nocolorbar
            #    colorbar();
            #end
            if nocolorbar
                plot!(colorbar = :none)
            end

            # custom dynamic display range?
            if cLim
                #set(gca, 'cLim', [cmin, cmax]);
                plot!(cLims = (cmin, cmax))
            end

            # symmetric?
            if symmetric
                if cLim
                    cmax = maximum([abs(cmin), abs(cmax)]);
                else
                    cmax = maximum(abs.(G.f[:]));
                end
                if cmax < 10*eps()
                    cmax = 1;
                end
                #set(gca, 'cLim', [-cmax, cmax]);
                plot!(cLims = (cmin, cmax))
            end
        end

        # title?
        if !isempty(G.name) && !notitle
            title!(G.name);
        end

        # verbatim?
        #if verbatim
        #    eval(verbtext);
        #end

    end

#set(gcf,'color','w');

#set(gca, 'FontSize',10);


end

"""
    ExtractPlotSpecs(G::GraphSig)

Extract the important information contained in the 'plotspecs' of a GraphSig object

### Input Argument
* `G::GraphSig`: an input GraphSig object
### Output Argument
*   symmetric       symmetrize the colorbar
*   markercolor     markercolor scheme
*   markerstrokealpha   capacity of marker stroke
*   markervaluevaries   if the marker color depends on the signal value
*   notitle         display a title
*   nocolorbar      display a colorbar
*   stemplot        use a stem plot
*   cLim            specify the dynamic display range
*   cmin            the min of the dynamic display range
*   cmax            the max of the dynamic display range
*   ptsize          the size of the nodes
*   linewide        the width of the lines in gplot
*   linecolor       the color of the lines (1D) / graph edges (2D & 3D)
*   verbatim        use other special instructions
*   verbtext        the specified special instructions
*   sortnodes       plot the signal values from smallest to largest in magnitude

"""
function ExtractPlotSpecs(G::GraphSig)
    plotspecs = G.plotspecs

    # symmetrize the colorbar?
    symmetric = occursin("symm", plotspecs)

    # marker color depends on the signal value?
    markervaluevaries = ~occursin("constant", plotspecs)

    # display a title?
    notitle  = occursin("notitle", plotspecs)

    # display a colorbar?
    nocolorbar = occursin("nocolorbar", plotspecs)

    # use a stem plot?
    stemplot = occursin("stem", plotspecs)

    # what should the dynamic display range be?
    temp, cLim = find_between(plotspecs,"cLim[","]")
    if cLim
        temp = split(temp,",")
        cmin = parse(Float64,temp[1])
        cmax = parse(Float64,temp[2])
    else
        temp, cLim = find_between(plotspecs,"cLim([","])")
        if cLim
            temp = split(temp,",")
            cmin = parse(Float64,temp[1])
            cmax = parse(Float64,temp[2])
        else
            cmin = []
            cmax = []
        end
    end

    # how big should the points be?
    if G.dim == 1
        ptsize = [2]
    else
        ptsize = [5]
    end

    temp, TF = find_between(plotspecs, "size[","]")
    if TF
        ptsize = parse(Float64,temp)
    end

    # how wide should the lines be
    if G.dim == 1
        linewide = 1
    else
        linewide = []
    end

    temp, TF = find_between(plotspecs,"width[","]")
    if TF
        linewide = parse(Float64,temp)
    end

    # what color should the lines (1D) / graph edges (2D & 3D) be ?
    if G.dim == 1
        linecolor = :blue
    else
        linecolor = :black
    end

    temp, TF = find_between(plotspecs,"color[","]")
    if TF
        linecolor = temp
    end


    # other plot instructions ?
    verbtext, verbatim = find_between(plotspecs, "verbatime{{","}}")

    # sort the nodes?
    sortnodes = false
    if occursin("sortnodes",plotspecs)
        sortnodes = true
    end

    #
    temp, TF = find_between(plotspecs,"markerstrokealpha[","]")
    if TF
        markerstrokealpha = parse(Float64,temp)
    else
        markerstrokealpha = 1.0
    end

    #
    temp, TF = find_between(plotspecs,"colorm[","]")
    if TF
        markercolor = Symbol(temp)
    else
        markercolor = :balance
    end


    return symmetric,markercolor,markerstrokealpha,markervaluevaries,notitle,nocolorbar,stemplot,cLim,cmin,cmax,ptsize,linewide,linecolor,verbatim,verbtext,sortnodes
end

function find_between(plotspecs::String,left::String,right::String)
# Given a string 'plotspecs', return str, where <left>str<right>
    TF = false;
    str = "";
    a = findfirst(left, plotspecs)
    if a!=nothing && length(plotspecs) > a[1]+length(left)-1
        b = findfirst(right, plotspecs[a[1]+length(left):end])[1]
        if b!=nothing
            str = plotspecs[a[1]+length(left):a[1]+length(left)+b[1]-2];
            TF = true;
        end
    end
    return str, TF
end
