# TWIG ReadMe

## Workflow

Carrying out a TWIG analysis is currently a two-step process. 

1. Run a Julia script to generate the Fisher Information Matrix eigenanalysis results for your system.
2. Visualize these results in R.

The Julia script requires the ParametricModels package, which you can acquire from the [BYU Modeling Hub](https://git.physics.byu.edu/Modeling/ParametricModels.jl).

## Example

Imagine you wish to analyze the subcritical pitchfork bifurcation definied by

$\dot{y} = r y - y^3 + \alpha_1 y^4 + \alpha_2 y^5 + ...$

which we know has a bifurcation when all parameters are zero, or $r=\alpha_i=0$. 

This system is analyzed with the `bifurc_pitch.jl` script, and the results are visualized in R using the `twig.r` script.

First, set the parameter values in a Julia struct

```
@parameterspace mutable struct Bif
   y0     = 0.5, identity
   r      = 0, identity
   alpha1 = 0, identity
   alpha2 = 0, identity
   alpha3 = 0, identity
   alpha4 = 0, identity
   alpha5 = 0, identity
end
```
Then, define the right-hand side of your equation. Note that we use y[1], since the ODE system can involve many equations if desired. The parameter set `ps` is defined as a variable of time `Bif`, the structure for holding all the parameters above.
```
function rhs(ps::Bif{T}, t, y, dy) where T<:Real
	dy[1]=ps.r*y[1]-y[1]^3+ps.alpha1*y[1]^4+ps.alpha2*y[1]^5+ps.alpha3*y[1]^6+ps.alpha4*y[1]^7+ps.alpha5*y[1]^8	
    nothing
end
```

Finally, set initial conditions, and run the model. It will generate several output files, including
- eigvec_pitchf.txt: a concatination of the all the eigenvector matrices evaluated at each value of $t_{max}$
- eigval_pitchf.txt: a matrix of all the eigenvalues associated with each eigenvector
- tmaxs.txt: a list of the $t_{max}$ values where the FIM was calculated and analyzed
- lsing_pitchf.txt: a concatination of all the left-singular value matrices

These results can then be analyzed in R, using `twig.R`. Simply copy-and-paste the entire `twig()` function into the console. Then run the following command:
```
> twig("~/bifurcations/","pitchf", pars=expression(y[0]==.5,r==0,alpha[1]==0,alpha[2]==0,alpha[3]==0,alpha[4]==0,alpha[5]==0), npart=2, formu=expression(dot(y)==r*y-y^3+alpha[1]*y^4+alpha[2]*y^5+alpha[3]*y^6+alpha[4]*y^7+alpha[5]*y^8))
```

This will generate the rainbow plot seen in the paper, and in `pitchf.png`.

## Further references

Good luck! For more information: 
- see our paper, available on the arXiv: https://arxiv.org/abs/2201.08301
- Or contact developer Christian N. K. Anderson by email: bifurcate@byu.edu
