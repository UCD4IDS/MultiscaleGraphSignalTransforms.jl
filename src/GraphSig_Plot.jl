using Plots
"""
    GraphSig_Plot(G::GraphSig; symmetric::Bool = false,
      markersize::Float64 = 2.,
      markercolor::Symbol = :balance,
      markershape::Symbol = :circle,
      markerstrokewidth::Float64 = 1.0,
      markerstrokealpha::Float64 = 1.0,
      markervaluevaries::Bool = true,
      linewidth::Float64 = 1.,
      linecolor::Symbol = :blue,
      linestyle::Symbol = :solid,
      clim::Tuple{Float64,Float64} = (0., 0.),
      notitle::Bool = false, nocolorbar::Bool = false, nolegend::Bool = true,
      stemplot::Bool = false, sortnodes::Bool = false)

Display a plot of the data in a GraphSig object

### Input Argument
*   `G::GraphSig`:        an input GraphSig object
*   `symmetric`           symmetrize the colorbar
*   `markersize`          the size of the nodes
*   `markercolor`         markercolor scheme
*   `markershape`         shape of marker
*   `markerstrokewidth`   width of marker stroke
*   `markerstrokealpha`   capacity of marker stroke
*   `markervaluevaries`   if the marker color depends on the signal value
*   `linewidth`           the width of the lines in gplot
*   `linecolor`           the color of the lines (1D) / graph edges (2D & 3D)
*   `linestyle`:          the style of line
*   `notitle`             display a title
*   `nocolorbar`          display a colorbar
*   `nolegend`            display a legend
*   `stemplot`            use a stem plot
*   `clim`                specify the dynamic display range
*   `sortnodes`           plot the signal values from smallest to largest in magnitude
"""
function GraphSig_Plot(G::GraphSig; symmetric::Bool = false,
      markersize::Any = 2.,
      markercolor::Symbol = :balance,
      markershape::Symbol = :circle,
      markerstrokewidth::Float64 = 1.0,
      markerstrokealpha::Float64 = 1.0,
      markervaluevaries::Bool = true,
      linewidth::Float64 = 1.,
      linecolor::Symbol = :blue,
      linestyle::Symbol = :solid,
      clim::Tuple{Float64,Float64} = (0., 0.),
      notitle::Bool = false, nocolorbar::Bool = false, nolegend::Bool = true,
      stemplot::Bool = false, sortnodes::Symbol = :normal)

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
    # set(0, 'DefaultFigureVisible', 'on');

    #################################
    ############# 1-D ###############
    #################################
    if G.dim == 1
        if G.length < 65 || stemplot
            # stem(G.xy, G.f, '-o','LineWidth',2,'Color',linecolor,'MarkerFaceColor',linecolor,'MarkerSize',markersize(1));
            # stem(G.xy, G.f, linetype = [:sticks :scatter], linewidth = 2, linecolor = linecolor, markercolor = linecolor, markersize = markersize)
            stem(G.f; x = G.xy, lw = linewidth, c = linecolor, ms = markersize, shape = markershape)
        else
            # plot(G.xy, G.f, '-','Color',linecolor,'LineWidth',linewidth);
            plot(G.xy, G.f, linetype = :line, linecolor = linecolor, linewidth = linewidth);
        end

        xmin = minimum(G.xy)
        xmax = maximum(G.xy)
        dx = xmax - xmin
        if dx < 10 * eps()
            dx = 1
        end
        xmin = xmin - 0.1 * dx
        xmax = xmax + 0.1 * dx

        ymin = minimum(G.f)
        ymax = maximum(G.f)
        dy = ymax - ymin
        if dy < 10 * eps()
            dy = 1
        end
        ymin = ymin - 0.1*dy
        ymax = ymax + 0.1*dy

        # symmetric?
        if symmetric
            ymax = 1.1 * max(ymax, -ymin)
            ymin = -ymax
        end

        xlims!(xmin, xmax)
        ylims!(ymin, ymax)

        # title?
        if ~isempty(G.name) && ~notitle
            title!(G.name)
        end

        # legend?
        if nolegend
            plot!(leg = false)
        end


        #################################
        ######### 2-D  or 3-D ###########
        #################################

    else

        #fig = figure('visible','on');

        # plot the graph
        gplot(G.W, G.xy; style = linestyle, color = linecolor, width = linewidth)
        #hold on

        # plot the nodes
        scatter_gplot!(G.xy;
                       marker = G.f, ms = markersize,
                       shape = markershape, mswidth = markerstrokewidth,
                       msalpha = markerstrokealpha, plotOrder = sortnodes,
                       c = markercolor)

        # colorbar?
        if nocolorbar
            plot!(colorbar = :none)
        end

        # custom dynamic display range?
        if clim != (0., 0.)
            plot!(clim = clim)
        end

        # symmetric?
        if symmetric
            if clim != (0., 0.)
                cmax = maximum(abs.(clim));
            else
                cmax = maximum(abs.(G.f[:]));
            end
            if cmax < 10*eps()
                cmax = 1;
            end
            plot!(cLims = (cmax, cmax))
        end


        # title?
        if !isempty(G.name) && !notitle
            title!(G.name);
        end

        # legend?
        if nolegend
            plot!(leg = false)
        end
    end
end
