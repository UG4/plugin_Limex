--------------------------------------------------------------------------------
--	tut08_nonlinear_conv_diff_using_self_coupling.lua
--
--	This tutorial is used to show tha ability to compute non-linear problems.
--	The problem itself will be a simple convection-diffusion equation that is
--	non-linear since the coefficients are non-linear.
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

local dim = util.GetParamNumber("-dim", 2) -- default dimension is 2.
InitUG(dim, AlgebraType("CPU", 1))
local gridName = util.GetParam("-grid", "richardson_radu.ugx")


local outFileNamePrefix = util.GetParam("-o", "distributed_domain_")
local numTimeSteps = util.GetParamNumber("--numTimeSteps", 10) -- default dimension is 10
local numPreRefs = util.GetParamNumber("--numPreRefs", 1)
local numTotalRefs = util.GetParamNumber("--numTotalRefs", 3)
local baseLevel= util.GetParamNumber("--baseLevel", 3)

-- Calculate the number of post-refs and make sure that the result makes sense.
local numPostRefs = numTotalRefs - numPreRefs
if numPostRefs < 0 then
	print("WARNING:\tnumPreRefs exceeds the number of total refinements.")
	print("\t\t\tNo refinement will be preformed after distribution.")
	numPostRefs = 0
end

util.CheckAndPrintHelp("Time-dependent non-linear convection diffusion example");


print("Loading domain from " .. gridName)
local dom = Domain()
LoadDomain(dom, gridName)


-- refinement
local refiner = GlobalDomainRefiner(dom)
for i = 1, numPreRefs do
	refiner:refine()
end

if util.DistributeDomain(dom) == false then
	print("Error while Distributing Domain. Aborting.")
	exit()
end
-- perform post-refinement
for i = 1, numPostRefs do
	refiner:refine()
end

-- Lets save the domain on each process
outFileName = outFileNamePrefix .. ProcRank() .. ".ugx"
SaveDomain(dom, outFileName)
print("Saved domain to " .. outFileName)



-- Check the subset handler if all subsets are given
sh = dom:subset_handler()
if sh:get_subset_index("Inner") == -1 then
	print("Domain does not contain subset 'Inner'. Aborting.")
	exit()
end

if sh:get_subset_index("Dirichlet1") == -1 then
	print("Domain does not contain subset 'Dirichlet1'. Aborting.")
	exit()
end

if sh:get_subset_index("Dirichlet2") == -1 then
  print("Domain does not contain subset 'Dirichlet2'. Aborting.")
  exit()
end




----------------------------------------------------
----------------------------------------------------
-- User function and User Data
----------------------------------------------------
----------------------------------------------------

-- Next we need a lot of User Functions and User Data to specify the 
-- convection-diffusion problem's coefficient. The data can be
-- set from the lua script.

local SiltLoamParams = {
  PhiS = 0.396,
  PhiR = 0.131,
  alpha = 0.423,
  n = 2.06,
  
  Ks= 4.96e-2,
  
  dt = 1.0/48.0,
  tEnd = 3.0/16.0,
  tDrain = 1.0/16.0,
}

local BeitNetofaClayParams = {
  PhiS = 0.446,
  PhiR = 0.0,
  alpha = 0.152,
  n = 1.17,
  
  Ks= 8.2e-4,
  
  dt = 1.0/3.0,
  tEnd = 3.0,
  tDrain = 1.0,
}


-- VanGenuchtenParams = SiltLoamParams
VanGenuchtenParams = BeitNetofaClayParams
-- defining \theta
function MySaturation(p)
  local n = VanGenuchtenParams.n
  local m = 1.0 - (1.0/n)
 -- print("p="..p)
 -- print("n="..n)
  local apn = math.pow(VanGenuchtenParams.alpha*p, n)
 -- print("apn="..apn)
  return math.pow(1.0+apn, (1.0/n)-1.0)
end

function Deriv_MySaturation(p)
  local n=VanGenuchtenParams.n
  local alpha=VanGenuchtenParams.alpha
  local apn = (alpha*p)^n
  local p1 = math.pow(1.0+apn, (1.0/n)-2.0)
  local p2 = math.pow(alpha*p, n-1.0)
  
  return alpha*(1.0-n)*p1*p2
end

-- defining \Phi
function MyPorosity(p)
  if (p<0) then return VanGenuchtenParams.PhiS end
  
  local thetaP=MySaturation(p)
 -- print("Saturation="..thetaP)
  return VanGenuchtenParams.PhiR + (VanGenuchtenParams.PhiS-VanGenuchtenParams.PhiR)*thetaP
end


-- defining derivatives of \Phi(p)
function Deriv_MyPorosity(p)

  if (p<0) then return 0.0 end
  
  local dtheta_dP=Deriv_MySaturation(p)
  return (VanGenuchtenParams.PhiS-VanGenuchtenParams.PhiR)*dtheta_dP
end

----------------------------------------------------
-- Conductivity & Diffusion Tensor
----------------------------------------------------

-- defining K(p)
function MyConductivity(p)
  local Ks= VanGenuchtenParams.Ks
  local n = VanGenuchtenParams.n
  
  if (p<0) then return Ks end

  local thetaP = MyPorosity(p)
  local thetaPow = math.pow(thetaP, n/(n-1.0)) 
  local val = math.pow(1.0-thetaPow, (n-1.0)/n)
  --print(thetaP)
  --print(val)
  local K = Ks*math.sqrt(thetaP)*(1.0-val)*(1.0-val)
 -- print("Conductivity="..K)
  return -K
end


-- defining Deriv K(p)
function Deriv_MyConductivity(p)

  if (p<0) then return 0 end
  
  local deltap = math.sqrt(math.abs(p)*1e-16+1e-16)
  return (MyConductivity(p+deltap)-MyConductivity(p))/deltap
end


function MyConductivityTensor(p)
	local K =  MyConductivity(p)
	return K, 0, 0, K
end

function Deriv_MyConductivityTensor(p)
	local DK = Deriv_MyConductivity(p)
	return 	DK, 0, 0, DK
end




----------------------------------------------------
-- Convection - Diffusion Element Discretization
----------------------------------------------------


print("Creating ApproximationSpace")
local approxSpace = ApproximationSpace(dom) -- creates new object
approxSpace:add_fct("p", "Lagrange", 1)     -- adds one function
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

local elemDisc = ConvectionDiffusion("p", "Inner", "fv1")

local LuaPorosity = LuaUserFunctionNumber("MyPorosity", 1);
LuaPorosity:set_deriv(0, "Deriv_MyPorosity");
LuaPorosity:set_input(0, elemDisc:value())

local LuaConductivity = LuaUserFunctionNumber("MyConductivity", 1);
LuaConductivity:set_deriv(0, "Deriv_MyConductivity");
LuaConductivity:set_input(0, elemDisc:value())

local MyGravity = ConstUserVector(0.0)
MyGravity:set_entry(dim-1, 1)

local MyFlux = ScaleAddLinkerVector()
MyFlux:add(LuaConductivity, elemDisc:gradient())
MyFlux:add(LuaConductivity, MyGravity)

--elemDisc:set_diffusion(LuaConductivity)	-- set linker diffusion matrix
elemDisc:set_mass(LuaPorosity)
elemDisc:set_flux(MyFlux)
--elemDisc:set_diffusion(1.0)
--elemDisc:set_velocity(ConstUserVector(1.0))
--elemDisc:set_mass_scale(1.0)
----------------------------------------------------
-- Dirichlet Boundary
----------------------------------------------------

function MyTrenchDirichletBnd2d(x, y, t)
  local tDrain = VanGenuchtenParams.tDrain
  if (t <= tDrain) then 
    return true, -2.0+2.2*(t/tDrain)
  end
  
	return true, 0.2
end

function MyConstDirichletBnd2d(x, y, t) 
return true, 1.0-y end

function MyInitialValues2d(x, y, t, si) 
return 1.0-y end

local dirichletBnd = DirichletBoundary()
dirichletBnd:add("MyTrenchDirichletBnd" .. dim .. "d", "p", "Dirichlet1")
dirichletBnd:add("MyConstDirichletBnd" .. dim .. "d", "p", "Dirichlet2")

----------------------------------------------------
-- Adding all Discretizations
----------------------------------------------------

-- Finally we create the discretization object which combines all the
-- separate discretizations into one domain discretization.
local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBnd)


----------------------------------------------------
----------------------------------------------------
-- Solver setup
----------------------------------------------------
----------------------------------------------------

-- we need a linear solver that solves the linearized problem inside of the
-- newton solver iteration. We use a geometric multi-grid method with
-- Jacobi smoothing and an LU base-solver.
local baseSolver = LU()
local ilut = ILUT()

local gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_solver(baseSolver)
gmg:set_base_level(baseLevel)
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_smoother(Jacobi(0.66))
gmg:set_cycle_type(1)
gmg:set_num_presmooth(3)
gmg:set_rap(true)
gmg:set_num_postsmooth(3)

local linSolver = LinearSolver()
linSolver:set_preconditioner(gmg)
linSolver:set_convergence_check(ConvCheck(100, 1e-12, 1e-12))
linSolver:set_compute_fresh_defect_when_finished(true)


-- Next we need a convergence check, that computes the defect within each 
-- newton step and stops the iteration when a specified creterion is fullfilled.
-- For our purpose the ConvCheck is sufficient. Please note,
-- that this class derives from a general IConvergenceCheck-Interface and
-- also more specialized or self-coded convergence checks could be used.
local newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

-- Within each newton step a line search can be applied. In order to do so an
-- implementation of the ILineSearch-Interface can be passed to the newton
-- solver. Here again we use the standard implementation.
local newtonLineSearch = StandardLineSearch()

-- Now we can set up the newton solver. We set the linear solver created above
-- as solver for the linearized problem and set the convergence check. If you
-- want to you can also set the line search.
local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)


-- since solver configurations can be quite complex, we print the configuration:
print("NewtonSolver configuration:")
print(newtonSolver:config_string())




----------------------------------------------------
----------------------------------------------------
-- Time loop
----------------------------------------------------
----------------------------------------------------


-- We create a grid function on the surface of our MultiGrid hierarchy.
u = GridFunction(approxSpace)


-- Next we have to initialize the solution with the start configuration.
Interpolate("MyInitialValues"..dim.."d", u, "p", 0.0);

-- In order to plot our time steps, we need a VTK writer. For time dependent
-- problems we start a time series. This is necessary to group the time 
-- series at the end of the time loop. We write the start solution at the beginning.
local out = VTKOutput()
out:print("Solution", u)



local concErrorEst = GridFunctionEstimator("p",2)
    
local endTime = VanGenuchtenParams.tEnd
local dt   = VanGenuchtenParams.dt    
local dtMin = 1e-2*dt
local dtMax = 1e-1*endTime



--[[ assemble matrix and rhs
local A = MatrixOperator()
local u = GridFunction(approxSpace)
local b = GridFunction(approxSpace)
domainDisc:assemble_linear(A, b)
SaveMatrixForConnectionViewer(u, A, "Stiffness.mat") 
--]]

if (false) then
---------------------------------------------
-- Standard (fixed size) Euler time stepping
---------------------------------------------

local time = 0.0
-- Now we create a time discretization. We use the theta-scheme. The time 
-- stepping scheme gets passed the domain discretization and will assemble
-- mass- and stiffness parts using the domain disc, then adding them in a
-- weighted way.
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit euler


-- Now we create an operator from the time discretization. We use the 
-- (nonlinear)-Operator interface.
local op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()
newtonSolver:init(op)
-- Since we use a one-step scheme, we need one extra solution vector to store
-- the old time step. This is done using the clone method. All previous 
-- solutions (in this case only one) are stored in the "PreviousSolution" 
-- object, that behaves like a queue of fixed size. We push the start solution
-- as the first old time step to our queue
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

numTimeSteps = math.ceil(endTime/dt)

-- Now we can start our time loop 
for step = 1, numTimeSteps do
  print("++++++ TIMESTEP " .. step .. "/"..numTimeSteps.." BEGIN ++++++")

  -- we choose the time step size to be constant here
  do_dt = dt
  
  -- setup time Disc for old solutions and timestep
  timeDisc:prepare_step(solTimeSeries, do_dt)
  
  -- prepare newton solver
  if newtonSolver:prepare(u) == false then 
    print ("Newton solver prepare failed at step "..step.."."); exit(); 
  end 
  
  -- apply newton solver
  if newtonSolver:apply(u) == false then
     print ("Newton solver apply failed at step "..step.."."); exit();
  end 

  -- compute the new (absolut) time
  time = solTimeSeries:time(0) + do_dt
  
  -- we write the newly computed time step to our time series
  out:print("Solution", u, step, time)
  
  -- get oldest solution
  oldestSol = solTimeSeries:oldest()

  -- copy values into oldest solution (we reuse the memory here)
  VecScaleAssign(oldestSol, 1.0, u)
  
  -- push oldest solutions with new values to front, oldest sol pointer is poped from end
  solTimeSeries:push_discard_oldest(oldestSol, time)

  print("++++++ TIMESTEP " .. step .. "  END ++++++");
end

-- non-linear adaptive time stepping 
--util.SolveLinearTimeProblem(u, domainDisc, linSolver, vtkFull, "AdaptiveSolution",
 --                   "ImplEuler", 1.0, 0.0, endTime, dt0, dtMin, dtMax);  
end


local adaptiveTimeStepConfig ={}
  
adaptiveTimeStepConfig["TOLERANCE"] = 0.01
adaptiveTimeStepConfig["REDUCTION"] = 0.5
adaptiveTimeStepConfig["ESTIMATOR"] = concErrorEst
      
adaptiveTimeStepConfig["DT"] = dt
adaptiveTimeStepConfig["DT_MIN"] = dtMin
adaptiveTimeStepConfig["DT_MAX"] = dtMax

if (false) then

-- non-linear adaptive time stepping 
util.SolveNonlinearProblemAdaptiveTimestep(u, domainDisc, newtonSolver, vtkFull, "AdaptiveSolution",
                    0.0, endTime, dt0, dtMin, dtMax, adaptiveTimeStepConfig);  
end


if (true) then
-------------------
-- LIMEX scheme
-------------------
 print (">> Setting up solver (3)")


local solverDesc = {
    name = "bicgstab",
    precond = {
      type = "gmg",
      approxSpace = approxSpace,
      smoother = {
        type ="ilu",
      },
      preSmooth = 2,      -- number presmoothing steps
      postSmooth = 2,    -- number postsmoothing steps
      baseLevel = 1,  -- 2 for 2d -- 1 for 3d
      rap=true
      
      
    },
    convCheck = {
      type ="standard",
      iterations = 100,
      absolute = 1e-12,
      reduction = 1e-9 }
  }



local concErrorEst = GridFunctionEstimator("p", 4)
concErrorEst:add("p", 2, 1, 1.0)



local adaptiveTimeStepConfig ={}

adaptiveTimeStepConfig["TOLERANCE"] = 0.01
adaptiveTimeStepConfig["REDUCTION"] = 0.5
adaptiveTimeStepConfig["ESTIMATOR"] = concErrorEst

adaptiveTimeStepConfig["DT"] = dt
adaptiveTimeStepConfig["DT_MIN"] = dt
adaptiveTimeStepConfig["DT_MAX"] = dt



local limexLSolver = {}
local limexNLSolver = {}

local limexConvCheck=ConvCheck(1, 5e-8, 1e-10, true)
limexConvCheck:set_supress_unsuccessful(true)

local lsolveCheck = ConvCheck(50, 1e-12, 1e-10)



local nstages = 3
for i=1,nstages do 

  limexLSolver[i] = util.solver.CreateSolver(solverDesc)
  limexLSolver[i]:set_convergence_check(lsolveCheck) 
    
  limexNLSolver[i] = NewtonSolver()
  limexNLSolver[i]:set_linear_solver(limexLSolver[i])
  limexNLSolver[i]:set_convergence_check(limexConvCheck)
  
  
  -- print(limexNLSolver[i])
end

-- print matrix (check for decoupling)
--limexNLSolver[1]:set_debug(dbgWriter)




local vtkobserver = VTKOutputObserver("LIMEx.vtk", out)

-- setup for time integrator
local limex = LimexTimeIntegrator(nstages)
for i=1,nstages do 
  limex:add_stage(i-1, i, domainDisc, limexNLSolver[i])
end

limex:add_error_estimator(concErrorEst)
limex:set_tolerance(0.01)
limex:set_time_step(dt)
limex:set_dt_min(dt*1e-3)
limex:set_dt_max(endTime*1e-1)
limex:set_increase_factor(2.0)
limex:set_stepsize_safety_factor(0.8)

limex:attach_observer(vtkobserver)

-- solve problem
print(">> Solve problem (3)")
limex:apply(u, endTime, u, 0.0)
out:write_time_pvd("ElderLimexMultistage.vtk", u)


end


print("")
print("done.")

