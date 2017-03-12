# 2D_HeatDistribution_MPI

Consider a square area that has known temperatures along each of its edges, and a single heat source at one boundary location. The temperature of the interior of an area will depend on the temperatures around it. A question that we would like to answer is: “How does the distribution of the heat in the area change depending on whether the heat source is near an edge or near the center of the area, and also depending on the temperature of the heat source?”
The approach to solving this problem is to divide the area into a fine mesh of points, h[i][j]. Temperatures at an inside point are taken to be the average of temperatures of four neighboring points.
