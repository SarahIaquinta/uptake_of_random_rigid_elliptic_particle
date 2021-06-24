# Cellular uptake of random elliptic particles

## Associated paper
The present code is the supplemental material associated to the paper [1]. 

## Dependencies
In order to make sure that you are able to run the code, please install the required versions of the libraries by executing the command bellow in your terminal.

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the requirements

```pip3 install -r requirements.txt```

## Figures
![beta](https://github.com/SarahIaquinta/uptake_of_random_rigid_elliptic_particle/blob/main/fig_beta.png)

## Tutorial
For given input values of the semi-major axis, semi-minor axis, adimensional lineic adhesion energy and adimensional membrane tension, this code allows you to:
- display the variation of the adimensional energy with respect to the wrapping degree
- determine the wrapping phase at equilibrium

To get the same results as the one presented in the paper [1], it is necessary to use the same input parameters, especially while setting the semi-major and semi-minor axes. Indeed, the perimeter of the particle should remain equal to 2*pi.

The values of f_list, sampling_points_membrane and sampling_points_circle should remain unchanged, as they result from convergence studies. 


To run the code in terminal, execute the following command to set the input parameters:

```sh
python3 uptake_of_random_rigid_elliptic_particle.py \
    --semi_major_axis 1 
    --semi_minor_axis 1 
    --sigma_bar 2
    --gamma_bar 10
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
A [GPL](https://tldrlegal.com/license/bsd-3-clause-license-(revised)) license is associated to this code, presented in the text file LICENSE.md.

## References
```
[1] @article{
        title={Cellular uptake of random rigid elliptic nanoparticles},
        author={Iaquinta S, Khazaie S, Fréour S, Jacquemin F, Blanquart C, Ishow E},
        journal={Physical Review Letters},
        year={2021}
        }
[2] @article{
        title={Ramanujan's Perimeter of an Ellipse},
        author={Villarino, Mark B},
        journal={arXiv preprint math/0506384},
        year={2005}
        }
```
