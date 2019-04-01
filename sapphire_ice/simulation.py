import firedrake as fe
import sapphire.simulations.convection_coupled_phasechange


def water_buoyancy(sim, temperature):
    """ Eq. (25) from @cite{danaila2014newton} """
    T = temperature
    
    T_anomaly_degC = fe.Constant(4.0293)
    
    rho_anomaly_SI = fe.Constant(999.972)
    
    w_degC = fe.Constant(9.2793e-6)
    
    q = fe.Constant(1.894816)
    
    M = sim.reference_temperature_range__degC
    
    def T_degC(T):
        """ T = T_degC/M """
        return M*T
    
    def rho_of_T_degC(T_degC):
        """ Eq. (24) from @cite{danaila2014newton} """
        return rho_anomaly_SI*(1. - w_degC*abs(T_degC - T_anomaly_degC)**q)
        
    def rho(T):
        
        return rho_of_T_degC(T_degC(T))
    
    beta = fe.Constant(6.91e-5)  # [K^-1]
    
    Gr = sim.grashof_number
    
    ghat = fe.Constant(-sapphire.simulation.unit_vectors(sim.mesh)[1])
    
    rho_0 = rho(T = 0.)
    
    return Gr/(beta*M)*(rho_0 - rho(T))/rho_0*ghat

    
def variational_form_residual(sim, solution):
    
    return sum(
    [r(sim = sim, solution = solution)
        for r in (
            sapphire.simulations.convection_coupled_phasechange.mass,
            lambda sim, solution: \
                sapphire.simulations.convection_coupled_phasechange.momentum(
                    sim = sim,
                    solution = solution,
                    buoyancy = water_buoyancy),
            sapphire.simulations.convection_coupled_phasechange.energy,
            sapphire.simulations.convection_coupled_phasechange.stabilization)])\
        *fe.dx(degree = sim.quadrature_degree)
    
    
class Simulation(sapphire.simulations.convection_coupled_phasechange.Simulation):

    def __init__(self,
            mesh,
            initial_values,
            dirichlet_boundary_conditions,
            quadrature_degree = 8,
            element_degree = 1,
            time_stencil_size = 2,
            output_directory_path = "output/"):
        
        self.reference_temperature_range__degC = fe.Constant(10.)
        
        super().__init__(
            variational_form_residual = variational_form_residual,
            mesh = mesh,
            initial_values = initial_values,
            dirichlet_boundary_conditions = dirichlet_boundary_conditions,
            quadrature_degree = quadrature_degree,
            element_degree = element_degree,
            time_stencil_size = time_stencil_size,
            output_directory_path = output_directory_path)
        
        self.stefan_number = self.stefan_number.assign(0.125)
        
        self.liquidus_temperature = self.liquidus_temperature.assign(0.)
        
        self.density_solid_to_liquid_ratio = \
            self.density_solid_to_liquid_ratio.assign(916.70/999.84)
        
        self.heat_capacity_solid_to_liquid_ratio = \
            self.heat_capacity_solid_to_liquid_ratio.assign(0.500)
        
        self.thermal_conductivity_solid_to_liquid_ratio = \
            self.thermal_conductivity_solid_to_liquid_ratio.assign(2.14/0.561)
        