import firedrake as fe 
import sapphire.test
import sapphire_ice.simulation


datadir = sapphire.test.datadir

def heat_driven_cavity_variational_form_residual(sim, solution):
    
    mass = sapphire.simulations.convection_coupled_phasechange.\
        mass(sim, solution)
    
    stabilization = sapphire.simulations.convection_coupled_phasechange.\
        stabilization(sim, solution)
    
    p, u, T = fe.split(solution)
    
    b = sapphire_ice.simulation.water_buoyancy(sim = sim, temperature = T)
    
    Pr = sim.prandtl_number
    
    _, psi_u, psi_T = fe.TestFunctions(sim.function_space)
    
    inner, dot, grad, div, sym = \
        fe.inner, fe.dot, fe.grad, fe.div, fe.sym
        
    momentum = dot(psi_u, grad(u)*u + b) \
        - div(psi_u)*p + 2.*inner(sym(grad(psi_u)), sym(grad(u)))
    
    energy = psi_T*dot(u, grad(T)) + dot(grad(psi_T), 1./Pr*grad(T))
    
    return mass + momentum + energy + stabilization
    
    
def dirichlet_boundary_conditions(sim):

    W = sim.function_space
    
    return [fe.DirichletBC(
        W.sub(1), (0., 0.), "on_boundary"),
        fe.DirichletBC(W.sub(2), sim.hot_wall_temperature, 1),
        fe.DirichletBC(W.sub(2), sim.cold_wall_temperature, 2)]
        
    
def initial_values(sim):
    
    print("Solving steady heat driven cavity to obtain initial values")
    
    Ra = 2.518084e6

    Pr = 6.99
    
    sim.reference_temperature_range__degC.assign(10.)
    
    sim.grashof_number = sim.grashof_number.assign(Ra/Pr)
    
    sim.prandtl_number = sim.prandtl_number.assign(Pr)
    
    w = fe.Function(sim.function_space)
    
    p, u, T = w.split()
    
    p.assign(0.)
    
    ihat, jhat = sim.unit_vectors()
    
    u.assign(0.*ihat + 0.*jhat)
    
    T.assign(sim.cold_wall_temperature)
    
    F = heat_driven_cavity_variational_form_residual(
        sim = sim,
        solution = w)*fe.dx(degree = sim.quadrature_degree)
    
    problem = fe.NonlinearVariationalProblem(
        F = F,
        u = w,
        bcs = dirichlet_boundary_conditions(sim),
        J = fe.derivative(F, w))
    
    solver = fe.NonlinearVariationalSolver(
        problem = problem,
        solver_parameters = {
                "snes_type": "newtonls",
                "snes_monitor": None,
                "ksp_type": "preonly", 
                "pc_type": "lu", 
                "mat_type": "aij",
                "pc_factor_mat_solver_type": "mumps"})
                
    def solve():
    
        solver.solve()
        
        return w
    
    w, _ = \
        sapphire.continuation.solve_with_bounded_regularization_sequence(
            solve = solve,
            solution = w,
            backup_solution = fe.Function(w),
            regularization_parameter = sim.grashof_number,
            initial_regularization_sequence = (
                0., sim.grashof_number.__float__()))
                
    return w
    

class WaterFreezingInCavitySimulation(
        sapphire_ice.simulation.Simulation):

    def __init__(self, *args, meshsize, **kwargs):
        
        self.reference_temperature_range__degC = fe.Constant(10.)
        
        self.hot_wall_temperature = fe.Constant(1.)
        
        self.cold_wall_temperature = fe.Constant(0.)
        
        super().__init__(
            *args,
            mesh = fe.UnitSquareMesh(meshsize, meshsize),
            initial_values = initial_values,
            dirichlet_boundary_conditions = dirichlet_boundary_conditions,
            **kwargs)
        
        self.stefan_number = self.stefan_number.assign(0.125)
        
        self.cold_wall_temperature = self.cold_wall_temperature.assign(-1.)
        
    
    
def freeze_water(endtime, s, tau, rx, nx, rt, nt, q, outdir = ""):
    
    mu_l__SI = 8.90e-4  # [Pa s]
    
    rho_l__SI = 999.84  # [kg / m^3]
    
    nu_l__SI = mu_l__SI/rho_l__SI  # [m^2 / s]
    
    t_f__SI = endtime  # [s]
    
    L__SI = 0.038  # [m]
    
    Tau = pow(L__SI, 2)/nu_l__SI
    
    t_f = t_f__SI/Tau
    
    """ For Kowalewski's water freezing experiment,
    at t_f__SI 2340 s, t_f = 1.44.
    """
    
    sim = WaterFreezingInCavitySimulation(
        quadrature_degree = q,
        element_degree = rx - 1,
        time_stencil_size = rt + 1,
        meshsize = nx,
        output_directory_path = str(outdir.join(
            "freeze_water/"
            + "s{0}_tau{1}/".format(s, tau)
            + "rx{0}_nx{1}_rt{2}_nt{3}/".format(rx, nx, rt, nt)
            + "q{0}/".format(q))))
    
    sim.timestep_size = sim.timestep_size.assign(t_f/float(nt))
    
    sim.solid_velocity_relaxation_factor = \
        sim.solid_velocity_relaxation_factor.assign(tau)
    
    sim.smoothing = sim.smoothing.assign(s)
    
    
    sim.solutions, _, = sim.run(endtime = t_f)
    
    
    print("Liquid area = {0}".format(sim.liquid_area))
    
    return sim
    
    
def test__validate__freeze_water__regression(datadir):
    
    sim = freeze_water(
        outdir = datadir,
        endtime = 2340.,
        s = 1./200.,
        tau = 1.e-12,
        rx = 2,
        nx = 24,
        rt = 2,
        nt = 4,
        q = 4)
    
    assert(abs(sim.liquid_area - 0.69) < 0.01)
    