# Cellular uptake of random elliptic particles

## Associated paper
The present code is the supplemental material associated to the paper [1]. 

## Dependencies
In order to make sure that you are able to run the code, please install the required versions of the libraries by executing the command bellow in your terminal.

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the requirements

```pip3 install -r requirements.txt```

## Figures

![coordinates](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/coordinates.png)
Figure 1: Definition of the coordinates system

![psi1](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/psi1_gif.gif)
Figure 2: Definition of psi 1 angle

![psi3](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/psi3_gif.gif)
Figure 3: Definition of psi 3 angle

![beta](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/beta_angles.png)
Figure 4: Definition of beta angles

![theta](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/theta_angles.png)
Figure 5: Definition of theta angles

![delta](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/figures/delta_angles.png)
Figure 5: Definition of delta angle

## Tutorial
For given input values of the aspect ratio and perimeter of the particle, adimensional lineic adhesion energy and adimensional membrane tension, this code allows you to:
- display the variation of the adimensional energy with respect to the wrapping degree
- determine the wrapping phase at equilibrium

The values of wrapping_list, sampling_points_membrane and sampling_points_circle are expected to be kept unchanged, as they result from convergence studies. 

To run the code in terminal, execute the command bellow to set the input parameters. 

```sh
python uptake_of_random_rigid_elliptic_particle.py \
    --r_bar 1 \
    --particle_perimeter 6.28 
```

To run it with the default values (r_bar = 1 and particle_perimeter = 2pi), execute only the following command:
```sh
python uptake_of_random_rigid_elliptic_particle.py
```

Remark: Depending on your python version, you might use a different command than "python", as "py", "py3" or "python3" for instance. 

It is also possible to run the code from any Python development environment. It will run the code written in the main section.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
A [GPL](https://tldrlegal.com/license/bsd-3-clause-license-(revised)) license is associated to this code, presented in the text file LICENSE.md.

## References
```
[1] @article{
        title={Influence of the mechanical and geometrical parameters on the cellular uptake of nanoparticles: a stochastic approach},
        author={Iaquinta S, Khazaie S, Fréour S, Jacquemin F, Blanquart C, Ishow E},
        journal={},
        year={2021}
        }
[2] @article{
        title={Ramanujan's Perimeter of an Ellipse},
        author={Villarino, Mark B},
        journal={arXiv preprint math/0506384},
        year={2005}
        }
```
