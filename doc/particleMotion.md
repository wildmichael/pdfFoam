# Particle Motion

The particle motion is modelled after [Jenny2001](references.md#Jenny2001). The
algorithm goes like follows:

1. Interpolate mean velocity to particle location:
   ![Mean velocity interpolation](http://quicklatex.com/cache3/09/ql_e8e555330a919ee76d2095c5cf272209_l3.png).
2. Perform a half-step

   ![First half step](http://quicklatex.com/cache3/dd/ql_1a2e688c6dd00d54569d56e7cc9abcdd_l3.png)

3. Interpolate mean velocity
   ![Mean velocity had half-step](http://quicklatex.com/cache3/33/ql_749931fc925bfaa56c794f3d1e04e533_l3.png),
   evaluate Langevin equation at new time level _n+1_ to estimate
   ![New particle velocity](http://quicklatex.com/cache3/ef/ql_a30441514b585d6a10871da4f8c229ef_l3.png)
   and perform mixing.
4. Perform full time step

   ![Full time step](http://quicklatex.com/cache3/af/ql_a72d66a8e58a3e5821e7d43ef13240af_l3.png)

To implement this algorithm efficiently im OpenFOAM, in particular to avoid
serious implementation issues at processor boundaries, the following modified
version of above algorithm is used:

1. Interpolate mean velocity to particle location:
   ![Mean velocity interpolation](http://quicklatex.com/cache3/09/ql_e8e555330a919ee76d2095c5cf272209_l3.png).
3. Perform a half-step

   ![First half step](http://quicklatex.com/cache3/4f/ql_0157094224861ef23c8754c53017424f_l3.png)

3. Interpolate mean velocity
   ![Mean velocity had half-step](http://quicklatex.com/cache3/33/ql_749931fc925bfaa56c794f3d1e04e533_l3.png),
   evaluate Langevin equation at new time level _n+1_ to estimate
   ![New particle velocity](http://quicklatex.com/cache3/ef/ql_a30441514b585d6a10871da4f8c229ef_l3.png)
   and perform mixing.
4. Compute the location of the particle at time level _n+1_
   according to the full time step

   ![Full time step](http://quicklatex.com/cache3/af/ql_a72d66a8e58a3e5821e7d43ef13240af_l3.png)

5. Calculate a &ldquo;tracking&rdquo; velocity moving the particle
   from ![Half step position](http://quicklatex.com/cache3/17/ql_d3fd2b2b6962885e09366d63301c2517_l3.png)
   to ![Full step position](http://quicklatex.com/cache3/ef/ql_9dace52dae425c44193b302b446ccaef_l3.png):

   ![Tracking velocity](http://quicklatex.com/cache3/f7/ql_8ffd41a9b178db99ad6c2214e88928f7_l3.png)

6. Perform the second half-step according to

   ![Second half step](http://quicklatex.com/cache3/b8/ql_e576cbb40644f2214c413eeaa745e6b8_l3.png)

The only issue I can currently think of is with strongly curved boundaries.
