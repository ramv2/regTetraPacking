<h2><u>Regular tetrahedron packing</u></h2>
An optimization scheme for packing regular tetrahedra. Implements the divide projection of the divide and concur algorithm discussed in detail here:
Y. Kallus, V. Elser, and S. Gravel, "Method for dense packing discovery", Physical Review E, vol. 82, no. 5, Nov. 2010.

Dependencies: 
1. org.apache.commons:commons-math3:3.6.1 (can be downloaded from Maven repository)
2. etomica.jar (see: https://github.com/etomica/etomica)
3. JUnit4 (can be downloaded from Maven repository)

Usage:
java -cp <path_to_regTetraPacking.jar>:<path_to_commons-math3-3.6.1.jar>:<path_to_etomica.jar> src.PackRegularTerahedra <int_value>

Here <...> denotes user input. For instance, the user must replace <path_to_commons-math3-3.6.1.jar> with the path to commons-math3-3.6.1.jar in her/his local machine. <int_value> denotes any input in the integer format.

<h2><u>Summary of approach</u></h2>
<ol>
<li> Obtain input 'N' denoting the number of regular tetrahedra to pack, from the user.</li>
<li> Initialize the possible length of edges to an array containing few values:
[0.01, 0.2, 0.1, 1, 2, 4, 8, 16, 100].</li>
<li> Loop through the array of edge lengths and let the current edge length
be l. Do the following:
<ul>
<li> Randomize the positions of N regular tetrahedra of length l in 3D
space.</li>
<li> Repeated pair-wise overlap detection and resolution based on
computing the divide projections for both tetrahedra. See
{@linkplain OverlapDetectionAndResolution}.</li>
<li> After convergence, loop through all position vectors to find the
range (i.e., max - min) of x-, y-, z- coordinates among all
tetrahedra vertices. This gives us the basis vectors and the volume
of the cell. </li>
<li> Evaluate the density and store the position vectors of all
tetrahedra as well as the basis vectors.</li>
</ul>
</li>
<li> Choose the l that results in the highest packing density as the result.</li>
</ol>
