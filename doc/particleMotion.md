# Particle Motion

The particle motion is modelled after [Jenny2001](references.md#Jenny2001). The
algorithm goes like follows:

1. Interpolate mean velocity to particle location:
   ![Mean velocity interpolation](https://latex.codecogs.com/svg.latex?\tilde{\vec{U}}\left(\vec{X}^{*^n}\right%29 )
2. Perform a half-step

   ![First half step](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+\frac{1}{2}}}=\vec{X}^{*^n}+\frac{\Delta&space;t}{2}\left(\tilde{\vec{U}}\left(\vec{X}^{*^n}\right%29+\vec{u}^{*^n}\right%29 )


3. Interpolate mean velocity
   ![Mean velocity at first half-step](https://latex.codecogs.com/svg.latex?\tilde{\vec{U}}\left(\vec{X}^{*^{n+\frac{1}{2}}}\right%29 )
   evaluate Langevin equation at new time level _n+1_ to estimate
   ![New particle velocity](https://latex.codecogs.com/svg.latex?\vec{u}^{*^{n+1}})
   and perform mixing.
4. Perform full time step

   ![Full time step](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+1}}=\vec{X}^{*^n}+\Delta&space;t\left(\tilde{\vec{U}}\left(\vec{X}^{*^{n+\frac{1}{2}}}\right%29+\frac{\vec{u}^{*^n}+\vec{u}^{*^{n+1}}}{2}\right%29 )

To implement this algorithm efficiently im OpenFOAM, in particular to avoid
serious implementation issues at processor boundaries, the following modified
version of above algorithm is used:

1. Interpolate mean velocity to particle location:
   ![Mean velocity interpolation](https://latex.codecogs.com/svg.latex?\tilde{\vec{U}}\left(\vec{X}^{*^n}\right%29 )
3. Perform a half-step

   ![First half step](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+\frac{1}{2}}}=\vec{X}^{*^n}+\frac{\Delta&space;t}{2}\left(\tilde{\vec{U}}\left(\vec{X}^{*^n}\right%29+\vec{u}^{*^n}\right%29 )

3. Interpolate mean velocity
   ![Mean velocity at first half-step](https://latex.codecogs.com/svg.latex?\tilde{\vec{U}}\left(\vec{X}^{*^{n+\frac{1}{2}}}\right%29 )
   evaluate Langevin equation at new time level _n+1_ to estimate
   ![New particle velocity](https://latex.codecogs.com/svg.latex?\vec{u}^{*^{n+1}})
   and perform mixing.
4. Compute the location of the particle at time level _n+1_
   according to the full time step

   ![Full time step](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+1}}=\vec{X}^{*^n}+\Delta&space;t\left(\tilde{\vec{U}}\left(\vec{X}^{*^{n+\frac{1}{2}}}\right%29+\frac{\vec{u}^{*^n}+\vec{u}^{*^{n+1}}}{2}\right%29 )

5. Calculate a &ldquo;tracking&rdquo; velocity moving the particle
   from ![Half step position](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+\frac{1}{2}}})
   to ![Full step position](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+1}}):

   ![Tracking velocity](https://latex.codecogs.com/svg.latex?\vec{U}^*_\text{track}=\frac{\vec{X}^{*^{n+1}}-\vec{X}^{*^{n+\frac{1}{2}}}}{\Delta&space;t/2})

6. Perform the second half-step according to

   ![Second half step](https://latex.codecogs.com/svg.latex?\vec{X}^{*^{n+1}}=\vec{X}^{*^{n+\frac{1}{2}}}+\frac{\Delta&space;t}{2}\vec{U}^*_\text{track})

The only issue I can currently think of is with strongly curved boundaries.
