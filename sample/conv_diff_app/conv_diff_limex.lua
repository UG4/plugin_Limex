--------------------------------------------------------------------------------
--[[!
-- \file sample/conv-diff-limex.lua
-- \ingroup app_convdiff
-- \{
-- \author Arne Naegel
-- \brief Illustrates limex usage based on conv_diff.lua by Andreas Vogel
-- \}
-- ]]
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")

local params = {}
local dim = util.GetParamNumber("-dim", 2, "dimension")
local numPreRefs = util.GetParamNumber("--numPreRefs", 0, "number of refinements before parallel distribution")
local numRefs    = util.GetParamNumber("--num-refs",    7, "number of refinements")
local startTime  = util.GetParamNumber("--start", 0.0, "start time")
local endTime    = util.GetParamNumber("--end", nil, "end time")
local dt         = util.GetParamNumber("--dt", 1e-2, "time step size")
local doVTK      = util.HasParamOption("--with-vtk")

params.tol     = util.GetParamNumber("--limex-tol", 1e-2, "time step size")
params.nstages = util.GetParamNumber("--limex-nstages", 2, "limex stages (2 default)")
params.limex_partial_mask = util.GetParamNumber("--limex-partial", 0, "limex partial (0 or 3)")
params.limex_debug_level = util.GetParamNumber("--limex-debug-level", 5, "limex debug level (integer)")



-- This scales the amount of diffusion of the problem
local eps       = util.GetParamNumber("--eps", 1e-0, "strength of diffusion")
local doVelocity = false

-- if withRedist == true, you may alter the redistribution behavior over the
-- parameters declared in util/load_balancing_util.lua.
-- Note that numPreRefs is ignored in this case.
local withRedist = util.HasParamOption("--withRedist", false, "Enables a more sophisticated distribution approach.")

GetLogAssistant():set_debug_level("LIB_LIMEX", params.limex_debug_level)


util.CheckAndPrintHelp("Time-dependent problem setup example\n(by Andreas Vogel)");

local gridName = nil
if dim == 2 then gridName = util.GetParam("-grid", "unit_square_01_quads_2x2.ugx")
else print("Dimension "..dim.." not supported."); exit(); end

local diffTime = 1.0*(1.0*1.0)/eps;
local convTime = 0.1 --math.pi/50.0;
endTime = endTime or math.min(convTime, diffTime);

print(" Selected Parameter:")
print("    numRefs      = " .. numRefs)
print("    numPreRefs   = " .. numPreRefs)
print("    startTime 	  = " .. startTime)
print("    endTime 		  = " .. endTime)
print("    dt 			    = " .. dt)
print("    eps          = " .. eps)
print("    grid         = " .. gridName)
print("    endTime      = " .. endTime)


print("    limex_nstages      = " .. params.nstages)
print("    limex_partial_mask = " .. params.limex_partial_mask)
print("    limex_tol          = " .. params.tol)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"Inner", "Boundary"}
local dom
if withRedist == true then
	dom = util.CreateDomain(gridName, 0, mandatorySubsets)
	balancer.RefineAndRebalanceDomain(dom, numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, mandatorySubsets)
end

print("\ndomain-info:")
print(dom:domain_info():to_string())

-- create Approximation Space
print(">> Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
--OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")


-- The coordinates (cx, cy) specify the rotation center of the cone
local cx = 0.5
local cy = 0.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
local ax = 0.25
local ay = 0.0

-- The parameter nu specifies the rotation velocity
local nu = 0.0-- 100

-- The parameter delta is a scaling factor influencing the steepness of the cone
delta = 1e-1  --1e-1
 
-- This is the exact solution for our problem
function exactSolution(x, y, t)
	local xRot = math.cos(nu*t) * (x-cx) - math.sin(nu*t) * (y-cy) 
	local yRot = math.sin(nu*t) * (x-cx) + math.cos(nu*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end
	
-- The velocity field
function Velocity(x, y, t)
	return	nu*(y - cx), nu*(cy - x)
end
	
-- The dirichlet condition
function DirichletValue(x, y, t, si)
	return true, exactSolution(x, y, t)
end


local vtk = nil 
local vtkObserver = nil
 
if (doVTK) then
   print("doVTK = true")
   vtk = VTKOutput();
   vtkObserver=VTKOutputObserver("MyFile.vtk", vtk)
end

-- post-processing (after each step)
function postProcess(u, step, time, currdt)
  local l2error = L2Error("exactSolution", u, "c", time, 4)
  local l2norm=L2Norm(u,"c", time)
  print("L2Error at t=\t"..time.."\t=\t"..l2error..", (rel=".. (l2error/l2norm).. ")\tL2Norm=\t"..l2norm)

  if (doVTK) then
     vtk:print("ConvDiffSol", u, step, time)
     local ref = u:clone()
     Interpolate("exactSolution", ref, "c", time)
     vtk:print("ConvDiffRef", ref, step, time)
  end
end



-- grid function debug writer
local dbgWriter = GridFunctionDebugWriter(approxSpace)
GetLogAssistant():set_debug_level("LIB_LIMEX", 4)
-- descriptor for linear solver
local solverDesc = {
    name = "bicgstab", -- "linear"
    precond = {
      type = "gmg",
      approxSpace = approxSpace,
      smoother = "sgs",
      rap=true,

      preSmooth=2,
      postSmooth=2,
      baseLevel = 2,
      baseSolver = "superlu",
      
--[[
      debugSolver = 
        {
          debug = true,
         type = "linear",    -- linear solver type ["bicgstab", "cg", "linear"]
            precond = { 
            type = "ilu",
             -- damping = 1.0,
          },
            
           
          convCheck = 
          {
            type    = "standard",
            iterations  = 10, -- number of iterations
            absolute  = 5e-12, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be larger than in newton section)
            reduction = 1e-8, -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be larger than in newton section)
            verbose   = true,   -- print convergence rates if true
          },
         },
--]]
      
    },
    convCheck = {
      type ="standard",
      maxSteps = 100,
      minDef = 1e-9,
      reduction = 1e-12 },
      
      
}



if (false) then
--------------------------------------------------------------------------------
--  Standard schemes
--------------------------------------------------------------------------------
print (">> Setting up Assembling (1)")
local elemDisc = ConvectionDiffusion("c", "Inner", "fv1")

elemDisc:set_upwind(FullUpwind())
elemDisc:set_diffusion(eps)
if (doVelocity) then
elemDisc:set_velocity("Velocity")
end
--elemDisc:set_partial_velocity(params.limex_partial_mask) -- 3 
elemDisc:set_partial_flux(0) 
elemDisc:set_partial_mass(0)  

local dirichletBND = DirichletBoundary()
dirichletBND:add("DirichletValue", "c", "Boundary")
--dirichletBND:add(0.0, "c", "Boundary")

local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)


print (">> Setting up solvers")
-- linear solver
util.solver.defaults.approxSpace  = approxSpace
local linSolver = util.solver.CreateSolver(solverDesc)

-- non-linear (Newton) Solver
local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(ConvCheck(10, 5e-8, 1e-10, true))
-- newtonSolver:set_debug(dbgWriter)

-- set initial value
print(">> Interpolation start values")
local u = GridFunction(approxSpace)
Interpolate("exactSolution", u, "c", startTime)


--[[ (non-linear) time loop
util.SolveNonlinearTimeProblem(u, domainDisc, newtonSolver, VTKOutput(), "Sol",
							   "ImplEuler", 1, startTime, endTime, dt); 
--]]

--[[ (linear) time loop
util.SolveLinearTimeProblem(u, domainDisc, linSolver, VTKOutput(), "Sol",
							"ImplEuler", 1, startTime, endTime, dt); 
--]]

--[[

LIMEX (outdated)

local dtMin = dt*1e-5;
local dtMax = endTime / 10.0;

local adaptiveTimeStepConfig ={
  -- ["ESTIMATOR"] =  GridFunctionEstimator("c", 2),
  ["ESTIMATOR"] =  Norm2Estimator(),
  ["TOLERANCE"] = params.tol,
  ["REDUCTION"] = 0.5,
  ["INCREASE"] = 4.0,
      
  ["DT"] = dt,
  ["DT_MIN"] = dtMin,
  ["DT_MAX"] = dtMax
}

local limexCheck = ConvCheck(1, 5e-6, 1e-10)
limexCheck:set_supress_unsuccessful(true)

newtonSolver:set_convergence_check(limexCheck)
-- newtonSolver:set_line_search(StandardLineSearch())

util.SolveNonlinearProblemAdaptiveLimex(u, domainDisc, newtonSolver, vtkFull, "Dummy",
                     startTime, endTime, dt, dtMin, dtMax, adaptiveTimeStepConfig, false, postProcess); 
--]]

else
print (">> Setting up Assembling (2)")
--------------------------------------------------------------------------------
--  Standard schemes
--------------------------------------------------------------------------------
local elemDisc ={}
local dirichletBND = {}
local domainDisc = {}


--local upwind = FullUpwind()
local upwind = NoUpwind()
-- setup for discretizations


if false then 
for i=1,nstages do 
  elemDisc[i] = ConvectionDiffusion("c", "Inner", "fv1")
  elemDisc[i]:set_upwind(upwind)
  elemDisc[i]:set_diffusion(eps)
  elemDisc[i]:set_velocity("Velocity")
  
  --elemDisc[i]:set_partial_velocity(params.limex_partial_mask) -- 3 
  --elemDisc[i]:set_partial_flux(0) 
  --elemDisc[i]:set_partial_mass(0)  
  
  dirichletBND[i] = DirichletBoundary()
  dirichletBND[i]:add("DirichletValue", "c", "Boundary")

  domainDisc[i] = DomainDiscretization(approxSpace)
  domainDisc[i]:add(elemDisc[i])
  domainDisc[i]:add(dirichletBND[i]) 
end 
else
   elemDisc = ConvectionDiffusion("c", "Inner", "fv1")
  elemDisc:set_upwind(upwind)
  elemDisc:set_diffusion(eps)
  elemDisc:set_velocity("Velocity")
  
  elemDisc:set_partial_velocity(params.limex_partial_mask) -- 3 
  elemDisc:set_partial_flux(0) 
  elemDisc:set_partial_mass(0)  
  
  dirichletBND = DirichletBoundary()
  dirichletBND:add("DirichletValue", "c", "Boundary")

  domainDisc = DomainDiscretization(approxSpace)
  domainDisc:add(elemDisc)
  domainDisc:add(dirichletBND) 
 end

 
local limexLSolver = nil
local limexNLSolver = nil

local limexConvCheck=ConvCheck(1, 5e-8, 1e-10, true)
limexConvCheck:set_supress_unsuccessful(true) 
 
if false then 
 limexLSolver = {}
 limexNLSolver = {}
for i=1,nstages do 

  limexLSolver[i] = util.solver.CreateSolver(solverDesc)
    
  limexNLSolver[i] = NewtonSolver()
  limexNLSolver[i]:set_linear_solver(limexLSolver[i])
  limexNLSolver[i]:set_convergence_check(limexConvCheck)
  
  print(limexNLSolver[i])
end
else
 limexLSolver = util.solver.CreateSolver(solverDesc)
  limexNLSolver = NewtonSolver()
 limexNLSolver:set_linear_solver(limexLSolver)
 limexNLSolver:set_convergence_check(limexConvCheck)
  print(limexNLSolver)
end


-- local refObserver = PlotRefOutputObserver("DirichletValue", vtk) -- now obsolete
local luaObserver = LuaCallbackObserver()

-- work-around (waiting for implementation of SmartPtr forward to lua...)
function luaPostProcess(step, time, currdt)
  print("LUAPostProcess: "..step..","..time..","..currdt)
  postProcess(luaObserver:get_current_solution(), step, time, currdt)
  return 0;
end
luaObserver:set_callback("luaPostProcess")

local tSteps = (endTime-startTime)/10.0
local dtlimex = tSteps 

local gridSize = 0.5*math.pow(0.5, numRefs)
if (doVelocity) then
  local tCFL = gridSize/50.0
  dtlimex = math.min(dtlimex, tCFL)
end
--  Euclidean (algebraic) norm
--local estimator = Norm2Estimator()  
--tol = 0.37/(gridSize)*tol


--print (estimator)
local limexEstimator = CompositeGridFunctionEstimator()
limexEstimator:add(L2ComponentSpace("c", 2))  

-- descriptor for integrator
local limexDesc = {

  nstages = params.nstages,
  steps = {1,2,3,4,5,6},
  domainDisc=domainDisc,
  nonlinSolver = limexNLSolver,
  
  tol = params.tol,
  dt = dtlimex,
  dtmin = 1e-9,
  
  rhoSafetyOPT = 0.25,
  
}


-- setup for time integrator
local limex = util.limex.CreateIntegrator(limexDesc)

limex:set_dt_min(1e-9)
limex:add_error_estimator(limexEstimator)
limex:set_increase_factor(2.0)

if (vtk) then 
   limex:attach_observer(vtkObserver)
end

limex:attach_observer(luaObserver)
--limex:attach_observer(refObserver)



limex:set_stepsize_greedy_order_factor(1.0)
limex:select_cost_strategy(LimexNonlinearCost())
-- limex:disable_matrix_cache()  -- recompute ()
limex:enable_matrix_cache() -- keep matrix 

print ("dtLimex   = "..dtlimex)
print ("hGrid     = "..gridSize)
print ("tolLimex  = "..params.tol)

dbgWriter:set_vtk_output(false)
dbgWriter:set_conn_viewer_output(true)
-- limexNLSolver:set_debug(dbgWriter)

limex:set_stepsize_greedy_order_factor(1.0)
limex:select_cost_strategy(LimexNonlinearCost())
-- limex:disable_matrix_cache()  -- recompute ()                                                                                                                                         
limex:enable_matrix_cache() -- keep matrix     

-- set initial value
print(">> Interpolating start values")
local u = GridFunction(approxSpace)
u:set(0.0)
Interpolate("exactSolution", u, "c", startTime)

-- solve problem

print(">> Peclet number:"..50.0*1.0/eps)
print(">> Grid Peclet number:"..50.0*gridSize/eps)
print(">> Solve problem")
local cstart=os.clock()
limex:apply(u, endTime, u, startTime)
local cend=os.clock()
print ("CDELTA=\t"..cend - cstart)

end

print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

FreeUserData()

-- end group app_convdiff
--[[!
\}
]]--
