"""
See LICENSE for details on how to use this code
"""
# Libraries
from functools import lru_cache
import argparse
from math import sin, cos, tan, atan, exp, pi, sqrt
from mpmath import csch, coth
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.signal


class MechanicalProperties:
    """
    A class to represent the mechanical properties of the cell membrane and the particle.

    Attributes:
        ----------
        gamma_bar: float
            adimensional linear adhesion energy between membrane and the particle
        sigma_bar: float
            adimensional membrane tension

    Methods:
        -------
        None

    """

    def __init__(self, gamma_bar, sigma_bar):
        """
        Constructs all the necessary attributes for the mechanical properties object.

        Parameters:
            ----------
            gamma_bar: float
                adimensional linear adhesion energy between membrane and the particle
            sigma_bar: float
                adimensional membrane tension

        Returns:
            -------
            None                    
        """
        self.gamma_bar = gamma_bar
        self.sigma_bar = sigma_bar


class ParticleGeometry:
    """
    A class to represent the particle.

    Attributes:
        ----------
        semi_minor_axis: float
            semi-minor axis of the elliptic particle
        semi_major_axis: float
            semi-major axis of the elliptic particle
        sampling_points_circle: int
            number of points to describe a circular particle with radius 1
        f: float
            wrapping degree 
        theta: float
            angle used as polar coordinate to define the ellipse, see figure * of readme

    Methods:
        -------
        define_particle_geometry_variables(self, f):
            Returns useful geometrical parameters to perform further calculations
        compute_r_coordinate(self, f, theta):
            Returns the value of the coordinate r given the position on the ellipse
        compute_z_coordinate(self, f, theta):
            Returns the value of the coordinate z given the position on the ellipse
        get_alpha_angle(self, f):
            Returns the value of the alpha angle (see figure *)
        compute_psi1_psi3_angles(self, f):
            Returns the curvature angles in the particle (see figure *)
        get_squared_dpsi_region3(self, f):
            Returns the values of dpsi3**2

    """
    def __init__(self, semi_minor_axis, semi_major_axis, sampling_points_circle):
        """
        Constructs all the necessary attributes for the particle object.

        Parameters:
            ----------
            semi_minor_axis: float
                semi-minor axis of the elliptic particle
            semi_major_axis: float
                semi-major axis of the elliptic particle
            sampling_points_circle: int
                number of points to describe a circular particle with radius semi_major axis.

        Returns:
            -------
            None                    
        """
        self.semi_minor_axis = semi_minor_axis
        self.semi_major_axis = semi_major_axis
        self.aspect_ratio = self.semi_major_axis / self.semi_minor_axis
        self.sampling_points_circle = sampling_points_circle
        h = ((self.semi_major_axis - self.semi_minor_axis) / \
             (self.semi_major_axis + self.semi_minor_axis))**2

        # computes the perimeter of the elliptic particle using Ramanujan's formula [2]
        self.particle_perimeter = pi * (self.semi_major_axis +  self.semi_minor_axis) * \
                                       (1 + 3*h / (10 + sqrt(4 - 3 * h)))
        self.effective_radius = self.particle_perimeter / (2 * pi)

        # amount of points to sample the particle.
        self.sampling_points_ellipse = self.sampling_points_circle * self.particle_perimeter / \
                                      (2*pi*self.semi_major_axis)

    @lru_cache(maxsize=10)                             
    def define_particle_geometry_variables(self, f):
        """
        Defines all the necessary variables to describe the elliptic particle
        for a given wrapping degree f.

        Parameters:
            ----------
            f: float
                wrapping degree (between 0 and 1)

        Returns:
            -------
            beta: float
                trigonometric angle at intersection between regions 1, 2r and 3 (see figure *)  
            beta_left: float
                trigonometric angle at intersection between regions 1, 2l and 3 (see figure *)       
            theta_list_region1: array
                trigonomic angle theta into the region 1 (see figure *)
            theta_list_region3: array
                trigonomic angle theta into the region 3 (see figure *)
            l1: float
                arclength of the region 1 (see figure *)
            l3: float
                arclength of the region 3 (see figure *)
            s_list_region1: array
                sampling of the arclength of the region 1 (see figure *)
            s_list_region3: array
                sampling of the arclength of the region 3 (see figure *)

        """
        beta = pi*f + 1.5 * pi
        beta_left = 3*pi - beta
        n3 = int(f * self.sampling_points_ellipse)
        n1 = int((1-f) * self.sampling_points_ellipse)
        theta_list_region1 = np.linspace(beta, 2*pi + beta_left, n1)
        theta_list_region3 = np.linspace(beta_left, beta, n3)
        l1 = self.particle_perimeter * (1 - f)
        l3 = self.particle_perimeter * f
        s_list_region1 = np.linspace(0, l1, n1)
        s_list_region3 = np.linspace(0, l3, n3)
        return beta, beta_left, theta_list_region1, theta_list_region3, s_list_region3, l3, s_list_region1, l1

    @lru_cache(maxsize=128)
    def compute_r_coordinate(self, f, theta):
        """
        Computes the value of the coordinate r given the position on the ellipse,
            depicted by the angle theta, for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree
            theta: float
                trigonomic angle in the ellipse

        Returns:
            -------
            r_coordinate: float
                r coordinate (see figure *)  
        """
        def compute_x_coordinate(t):
            """
            Computes the value of the coordinate x given the position on the ellipse,
                depicted by the angle theta

            Parameters:
                ----------
                t: float
                    angle in the ellipse

            Returns:
                -------
                x_coordinate: float
                    x coordinate (see figure *)  
            """
            x_coordinate = self.semi_major_axis * cos(t)
            return x_coordinate

        _, beta_left, _, _, _, _, _, _ = self.define_particle_geometry_variables(f)
        r_coordinate = compute_x_coordinate(theta) - compute_x_coordinate(beta_left)
        return r_coordinate

    @lru_cache(maxsize=128)
    def compute_z_coordinate(self, f, theta):
        """
        Computes the value of the coordinate z given the position on the ellipse,
            depicted by the angle theta, for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree
            theta: float
                trigonomic angle in the ellipse

        Returns:
            -------
            z_coordinate: float
                z coordinate (see figure *)  
        """
        def compute_y_coordinate(t):
            """
            Computes the value of the coordinate y given the position on the ellipse,
                depicted by the angle theta

            Parameters:
                ----------
                t: float
                    angle in the ellipse

            Returns:
                -------
                y_coordinate: float
                    y coordinate (see figure *)  
            """
            y_coordinate = particle.semi_minor_axis * sin(t)
            return y_coordinate

        _, beta_left, _, _, _, _, _, _ = self.define_particle_geometry_variables(f)
        z_coordinate = compute_y_coordinate(theta) - compute_y_coordinate(beta_left)
        return z_coordinate

    @lru_cache(maxsize=10)
    def get_alpha_angle(self, f):
        """
        Computes the value of the alpha angle (see figure *), for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree

        Returns:
            -------
            alpha: float
                curvature angle at intersection between regions 1, 2l and 3 (see figure *)  

        """
        psi_list_region1, _ = self.compute_psi1_psi3_angles(f)
        alpha = psi_list_region1[0]
        return alpha

    @lru_cache(maxsize=10)
    def compute_psi1_psi3_angles(self, f):
        """
        Computes the curvature angles in the particle (see figure *),
            for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree

        Returns:
            -------
            psi_list_region1: list
                psi angle in region 1 (see figure *)
            psi_list_region3: list
                psi angle in region 3 (see figure *)

        """
        beta, beta_left, theta_list_region1, theta_list_region3, s_list_region3, _, s_list_region1, _ \
            = self.define_particle_geometry_variables(f)
        x_bl = self.semi_major_axis * cos(beta_left)

        def compute_psi_from_r_z(theta):
            """
            Computes the curvature angle given the position in the ellipse,
                depicted by theta (see figure *), using r and z coordinates
            
            Parameters:
                ----------
                theta: float
                    trigonomic angle in the ellipse

            Returns:
                -------
                delta: float
                    angle between the tangent to the particle and horizontal (see figure *)

            """

            r_elli = self.compute_r_coordinate(f, theta)
            z_elli = self.compute_z_coordinate(f, theta)

            def compute_tangent_to_ellipse_at_rtan_and_theta(r_tan):
                """
                Computes the position of the tangent to the particle (z with respect to r)
                
                Parameters:
                    ----------
                    r_tan: float
                        r coordinates where the tangent equation is evaluated

                Returns:
                    -------
                    z_tan: float
                        z coordinate of the tangent at r = r_tan

                """

                r_elli = self.compute_r_coordinate(f, theta)
                z_elli = self.compute_z_coordinate(f, theta)

                # managing possible singularities
                if r_elli == (-x_bl - self.semi_major_axis):
                    r_elli = r_elli + 0.01 * self.semi_major_axis
                elif r_elli == (-x_bl + self.semi_major_axis):
                    r_elli = r_elli - 0.01 * particle.semi_major_axis
                slope1 = 1
                """
                Returns:
                    Eq of the tangent to the ellipse at theta = beta
                """
                # depending on the side of the particle, the tangent's slope is positive or negative
                if theta > pi:
                    slope1 = -1
                dz = -self.semi_minor_axis * (r_elli + x_bl) / (self.semi_major_axis **2) / \
                    sqrt(1 - ((r_elli + x_bl)/self.semi_major_axis)**2)
                z_tan = z_elli  + slope1 * dz * (r_tan - r_elli)
                return z_tan

            r_tan_list = np.linspace(0.5 * r_elli - 0.5, 1.5 * r_elli + 0.5, 10)
            r1 = max(r_tan_list)
            z1 = compute_tangent_to_ellipse_at_rtan_and_theta(r1)
            delta = atan(abs((z1 - z_elli)/(r1 - r_elli)))
            return delta

        delta_list_region1 = [compute_psi_from_r_z(t) for t in theta_list_region1]
        delta_list_region3 = [compute_psi_from_r_z(t) for t in theta_list_region3]
        psi_list_region1 = np.zeros_like(s_list_region1)
        psi_list_region3 = np.zeros_like(s_list_region3)

        if f < 0.5:  # psi angle is defined differently depending on the position on the particle
            for i in range(len(s_list_region3)):
                theta = theta_list_region3[i]
                delta = delta_list_region3[i]
                if theta < 1.5*pi:
                    psi = 2*pi - delta
                elif theta == 1.5*pi:
                    psi = 2*pi
                elif theta <= beta:
                    psi = 2*pi + delta
                psi_list_region3[i] = psi
            for i in range(len(s_list_region1)):
                theta = theta_list_region1[i]
                delta = delta_list_region1[i]
                if theta <= 2*pi:
                    psi = delta
                elif theta <= 2*pi + pi/2:
                    psi = pi - delta
                elif theta <= 3*pi:
                    psi = pi + delta
                else:
                    psi = 2*pi - delta
                psi_list_region1[i] = psi
        else:
            for i in range(len(s_list_region3)):
                theta = theta_list_region3[i]
                delta = delta_list_region3[i]
                if theta <pi:
                    psi = pi + delta
                elif theta == pi:
                    psi = 1.5*pi
                elif theta < 1.5*pi:
                    psi = 2*pi - delta
                elif theta == 1.5*pi:
                    psi = 2*pi        
                elif theta < 2*pi:
                    psi = 2*pi + delta
                elif theta == 2*pi:
                    psi = 2.5*pi
                elif theta <=beta:
                    psi = 3*pi - delta
                psi_list_region3[i] = psi

            for i in range(len(s_list_region1)):
                theta = theta_list_region1[i]
                delta = delta_list_region1[i]
                if theta < 2.5*pi:
                    psi = pi - delta
                elif theta == 2.5 * pi:
                    theta = pi
                elif theta <= 2*pi + beta_left:
                    psi = pi + delta
                psi_list_region1[i] = psi

        return psi_list_region1, psi_list_region3

    def get_squared_dpsi_region3(self, f):
        """
        Computes the values of dpsi3**2, necessary to evaluate the bending energy
            of the region 3, for a given wrapping degree f
        
        Parameters:
            ----------
            f: float
                wrapping degree

        Returns:
            -------
            squared_dpsi_list_region3: list
                dpsi angle power 2 in region 3 (see figure *)

        """
        _, _, _, _, s_list_region3, _, _, _ = self.define_particle_geometry_variables(f)
        _, psi_list_region3 = self.compute_psi1_psi3_angles(f)
        ds = s_list_region3[1] - s_list_region3[0]
        # computes the derivative of psi in region 3 using finite differences method. The value of
        # ds was set after a convergence study.
        dpsi3_list_region3 = [(psi_list_region3[i+1] - psi_list_region3[i])/ds 
                               for i in range(0, len(psi_list_region3)-1)]
        squared_dpsi_list_region3 = [p**2 for p in dpsi3_list_region3]
        return squared_dpsi_list_region3


class MembraneGeometry:
    """
    A class to represent the membrane object.

    Attributes:
        ----------
        particle: class
            ParticleGeometry class
        mechanics: class
            MechanicalProperties class
        sampling_points_membrane: int
            number of points to sample the regions 2r and 2l

    Methods:
        -------
        compute_r2r_r2l_z2r_z2l_from_analytic_expression(self, f, particle, mechanics):
            Returns r and z coordinates in the regions 2r and 2l

    """    
    def __init__(self, particle, sampling_points_membrane):
        """
        Constructs all the necessary attributes for the membrane object.

        Parameters:
            ----------
            particle: class
                ParticleGeometry class
            sampling_points_membrane: int
                number of points to sample the regions 2r and 2l

        Returns:
            -------
            None                    
        """
        self.sampling_points_membrane = sampling_points_membrane
        self.l2 = 20 * particle.effective_radius
        S2a = np.linspace(0, (self.l2 / 2), int((0.8*self.sampling_points_membrane) + 1))
        S2b = np.linspace(1.2 * self.l2 / 2, self.l2, int((0.2 * self.sampling_points_membrane)))
        self.S2 = np.concatenate((S2a, S2b), axis=None)

    @lru_cache(maxsize=10)
    def compute_r2r_r2l_z2r_z2l_from_analytic_expression(self, f, particle, mechanics):
        """
        Computes the r and z coordinates to describe the regions 2r and 2l,
            for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree
            particle: class
                ParticleGeometry class
            mechanics: class
                MechanicalProperties class
        Returns:
            -------
            r2r: list
                r coordinate in the region 2r
            z2r: list
                z coordinate in the region 2l
            r2l: list
                r coordinate in the region 2r
            z2l: list
                z coordinate in the region 2l
        """
        _, _, _, theta_list_region3, _, _, _, _ = particle.define_particle_geometry_variables(f)
        alpha = particle.get_alpha_angle(f)
        r2r_0 = particle.compute_r_coordinate(f, theta_list_region3[-1])
        z2r_0 = particle.compute_z_coordinate(f, theta_list_region3[-1])
        r2r = np.zeros_like(self.S2)
        z2r = np.zeros_like(self.S2)
        sigma = mechanics.sigma_bar

        for i in range(1, len(self.S2)):
            s = self.S2[i]
            r = r2r_0 + s - sqrt(2 / sigma) * (1 - cos(alpha)) / \
                (coth(s * sqrt(0.5*sigma)) + cos(alpha * 0.5))
            z = z2r_0 + sqrt(8 / sigma) * sin(0.5*alpha) * \
                (1 - (csch(s *sqrt(0.5*sigma))) / (coth(s* sqrt(0.5*sigma)) + cos(0.5*alpha)))
            r2r[i] = r
            z2r[i] = z

        r2r[0] = r2r_0
        z2r[0] = z2r_0
        r2l = np.array([r2r[0] - r2r[s] for s in range(len(self.S2))])
        z2l = z2r
        return r2r, z2r, r2l, z2l


def compute_adimensioned_energy_variation(f, particle, mechanics, membrane):
    """
    Computes the adimensional energy variation between a given wrapping degree f
        and the initial state at wero wrapping, as expressed in [1]

    Parameters:
        ----------
        f: float
            wrapping degree
        particle: class
            ParticleGeometry class
        mechanics: class
            MechanicalProperties class
        membrane: class
            MembraneGeometry class

    Returns:
        -------
        adimensional_total_energy_variation: float
            computed value of the adimensional total energy variation                     
    """
    _, _, _, _, s_list_region3, l3, _, _ = particle.define_particle_geometry_variables(f)
    r2r, _, r2l, _ = membrane.compute_r2r_r2l_z2r_z2l_from_analytic_expression(f,
                                                                               particle,
                                                                               mechanics)

    def bending_energy3():
        """
        Args:
            f: wrapping degree
            particle: Parameters class        
        Returns:
            float - the variation of the bending energy in the region 3 (bending of the membrane)
            between the state at wrapping f and the initial state: E(f) - E(0)
        """   
        dpsi3_list_region3 = particle.get_squared_dpsi_region3(f)
        energy = scipy.integrate.simps(dpsi3_list_region3, s_list_region3[0:-1])
        return energy

    def bending_energy2d():
        """
        Args:
            f: wrapping degree
            particle: Parameters class        
        Returns:
            float - the variation of the bending energy in the region 2d (free membrane)
            between the state at wrapping f and the initial state: E(f) - E(0)
        """     
        a = particle.get_alpha_angle(f)
        t = tan(0.25 * a)
        t2 = t**2
        b = sqrt(mechanics.sigma_bar / 2)
        E2 = - 8* b / particle.effective_radius * t2 * \
             ((1/(t2 + exp(2* b *membrane.l2 / particle.effective_radius))) - (1/(t2 + 1)))
        return E2

    Eb3 =  0.25 * particle.effective_radius * bending_energy3()
    Eb2d = 0.25 * particle.effective_radius * bending_energy2d()
    Eb2g = Eb2d
    Eb = Eb2d + Eb2g + Eb3
    Eadh = - mechanics.gamma_bar * l3 * 0.25 / particle.effective_radius
    Etens =  mechanics.sigma_bar * ( l3 + 2*membrane.l2 - (r2r[-1] - r2l[-1])) * \
             0.25 / particle.effective_radius
    adimensional_total_energy_variation = Eb + Eadh + Etens
    return adimensional_total_energy_variation


def plot_energy(particle, mechanics, membrane):
    """
    Plots the evolution of the adimensional variation of energy during wrapping

    Parameters:
        ----------
        particle: class
            ParticleGeometry class
        mechanics: class
            MechanicalProperties class
        membrane: class
            MembraneGeometry class

    Returns:
        -------
        None                 
    """
    energy_list = np.array([compute_adimensioned_energy_variation(f,
                                                                  particle,
                                                                  mechanics,
                                                                  membrane) for f in f_list])
    plt.figure()
    plt.plot(f_list,
             energy_list,
             '-k',
             label = '$\overline{r} = $' + str(np.round(particle.aspect_ratio, 2)) + 
                     ' ; $\overline{\gamma} = $' + str(mechanics.gamma_bar) + 
                     ' ; $\overline{\sigma} = $' + str(mechanics.sigma_bar))
    plt.xlabel('wrapping degree f [-]')
    plt.ylabel('$\Delta E [-]$')
    plt.title('$\Delta E(f)$')
    plt.legend()
    plt.xlim((0, 1))


def identify_wrapping_phase(particle, mechanics, membrane):
    """
    Identifies the wrapping phase following the process introduced in [1]

    Parameters:
        ----------
        particle: class
            ParticleGeometry object
        mechanics: class
            MechanicalProperties object
        membrane: class
            MembraneGeometry object

    Returns:
        -------
        wrapping_phase_number: float
            phase number (1, 2 or 3)
        wrapping_phase: str
            the wrapping phase as an intelligible string         
    """
    def get_local_energy_minima():
        """
        Identifies the local minima of energy variation with their location

        Parameters:
            ----------
            None

        Returns:
            -------
            energy_list:list
                list of adimensional energy variation during wrapping
            min_energy_list: list
                list of local minimna of adimensional energy with respect to wrapping degree f
            f_min_energy_list: list
                list of wrapping degree values where local minima of adimensional energy are reached
        """
        energy_list = np.array([compute_adimensioned_energy_variation(f,
                                                                      particle,
                                                                      mechanics,
                                                                      membrane) for f in f_list])
        min_energy_index_list = scipy.signal.argrelextrema(energy_list, np.less)

        # check if the minimum is reached for f_list[-1]
        if energy_list[-1] <  energy_list[-2]:
            min_energy_index_list = np.concatenate((min_energy_index_list,
                                                    np.array([-1])), axis=None)

        # check if the minimum is reached for f_list[0]
        if energy_list[0] < energy_list[1]:
            min_energy_index_list = np.concatenate((np.array(f_list[0]),
                                                    min_energy_index_list), axis=None)

        min_energy_list = [energy_list[int(k)] for k in min_energy_index_list]
        f_min_energy_list = [f_list[int(k)] for k in min_energy_index_list]

        return energy_list, min_energy_list, f_min_energy_list
    
    _, _, f_min_energy_list = get_local_energy_minima()
    f_eq = f_min_energy_list[0]
    wrapping_phase_number = 0
    wrapping_phase = '0'

    if f_eq < 0.2:  # check if wrapping phase is phase 1, according to [1]
        wrapping_phase_number = 1
        wrapping_phase = 'no wrapping'
    else: 
        r2r, _, r2l, _ = membrane.compute_r2r_r2l_z2r_z2l_from_analytic_expression(f_eq,
                                                                                   particle,
                                                                                   mechanics)
        intersection_membrane = min(r2r) - max(r2l)
        wrapping_phase_number = 3 if intersection_membrane < 0 else 2
        wrapping_phase = 'full wrapping' if intersection_membrane < 0 else "partial wrapping"
    
    def plot_geometry(f, particle, mechanics, membrane): 
        """
        Plots the morphology of the interaction between the cell and the particle,
            for a given wrapping degree f

        Parameters:
            ----------
            f: float
                wrapping degree
            particle: class
                ParticleGeometry class
            mechanics: class
                MechanicalProperties class
            membrane: class
                MembraneGeometry class

        Returns:
            -------
            None                    
        """      
        _, _, theta_list_region1, theta_list_region3, _, _, _, _ = \
        particle.define_particle_geometry_variables(f)
        r1 = [particle.compute_r_coordinate(f, theta) for theta in theta_list_region1]
        r3 = [particle.compute_r_coordinate(f, theta) for theta in theta_list_region3]
        z1 = [particle.compute_z_coordinate(f, theta) for theta in theta_list_region1]
        z3 = [particle.compute_z_coordinate(f, theta) for theta in theta_list_region3]
        r2r, z2r, r2l, z2l = \
        membrane.compute_r2r_r2l_z2r_z2l_from_analytic_expression(f, particle, mechanics)
        plt.figure()
        plt.plot(r1, z1, '-b', label = 'region 1')
        plt.plot(r2r, z2r, '-k', label = 'region 2')
        plt.plot(r2l, z2l, '-k')
        plt.plot(r3, z3, '-r', label = 'region 3')
        plt.xlabel('r (x100) (nm)')
        plt.ylabel('z (x100) (nm)')
        plt.ylim((-4, 4))
        plt.title('wrapping morphology')
        plt.gca().set_aspect('equal' , adjustable='box')
        plt.legend()    

    plot_geometry(f_eq, particle, mechanics, membrane)
    plt.title('wrapping morphology at equilibrium: ' + wrapping_phase)

    return wrapping_phase_number, wrapping_phase


def parse_arguments():
    """
    Parses arguments to run the code in terminal

    Parameters:
        ----------
        None

    Returns:
        -------
        args: class
            #TODO complete here
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gamma_bar', required=False, default=10.0, type=float,
        help='adimensional lineic adhesion between the membrane and the particle. Default value = 10.')
    parser.add_argument('-s', '--sigma_bar', required=False, default=2.0, type=float,
        help='adimensional membrane tension. Default value = 2.')
    parser.add_argument('-a', '--semi_major_axis', required=False, default=1.0, type=float,
        help='semi-major axis of the elliptic particle. Default value = 1.')
    parser.add_argument('-b', '--semi_minor_axis', required=False, default=1.0, type=float,
        help='semi-minor axis of the elliptic particle. Default value = 1.')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    f_list = np.arange(0.03, 0.97, 0.003125)

    particle = ParticleGeometry(semi_minor_axis=args.semi_minor_axis,
                                semi_major_axis=args.semi_major_axis,
                                sampling_points_circle=300)

    mechanics = MechanicalProperties(gamma_bar=args.gamma_bar,
                                     sigma_bar=args.sigma_bar)

    membrane = MembraneGeometry(particle,
                                sampling_points_membrane=100)

    plot_energy(particle, mechanics, membrane)
    wrapping_phase_number, wrapping_phase = identify_wrapping_phase(particle, mechanics, membrane)
    print('wrapping phase at equilibrium: ', wrapping_phase)
    plt.show()