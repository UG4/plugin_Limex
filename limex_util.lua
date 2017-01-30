-- Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
-- Author: Arne Naegel
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


--[[ (sample) solver descriptor

local limexDescSerial = {

  nstages = 2,
  steps = {1,2,3},
  domainDisc=domainDisc
  nonlinSolver = nlSolver
 
 --linSolver = linSolver,
  
  nthreads = 1
  tol = tol,
  dt = dtlimex,
  dtmin = 1e-9,
  
}


local limexDescParallel = {

  nstages = 2,
  steps = {1,2,3},
  domainDisc={domainDisc[1], domainDisc[2], domainDisc[3]}
  nonlinSolver = {nlSolver[1], nlSolver[2], nlSolver[3]}
 --lSolver = limexLSolver,
  
  nthreads = 2
  tol = tol,
  dt = dtlimex,
  dtmin = 1e-9,
  
}

--]]

ug_load_script("util/table_desc_util.lua")
ug_load_script("util/table_util.lua")
ug_load_script("util/solver_util.lua")

util = util or {}
util.limex = util.limex or {}

util.limex.defaultDesc = util.limex.defaultDesc or 
{
    nstages = 2,
    steps = {1,2,3,4,5,6,7,8,9,10},
    nthreads = 1,
    tol = 0.001,
    
    makeConsistent = false,
}

-- aux function creating a solver
function util.limex.CreateLimexSolver(nlsolverDesc, solverutil)
     local limexSolverDesc = {
        type = "newton",
        lineSearch = "none",
        convCheck = nlsolverDesc.convCheck,
     }
     
    return util.solver.CreateSolver(limexSolverDesc, solverutil)
  
end

-- aux function creating an integrator
function util.limex.CreateIntegrator(limexDesc)

-- max number of stages [scalar]
local nstages = limexDesc.nstages
if ((type(limexDesc.steps) ~= "table")) then print ("ERROR: Requires array of steps!") return nil
end

-- distribution of steps [scalar or string???]
local nsteps = table.getn(limexDesc.steps)
if ((nsteps < nstages)) then print ("ERROR: Array too short!") return nil 
end

print(limexDesc.nonlinSolver:config_string())

local ndiscs = 0 -- discretization(s) [object or table of objects]
if (type(limexDesc.domainDisc) == "table") then ndiscs = table.getn(limexDesc.domainDisc) 
end

local nsolvers = 0 -- solver(s) [object or table of objects]
if (type(limexDesc.nonlinSolver) == "table") then nsolvers = table.getn(limexDesc.nonlinSolver)
end


if (ndiscs ~= nsolvers) then 
  print ("Discs: "..ndiscs..","..nsolvers)
  print ("ERROR: domainDisc and nonlinSolver must match in type and number of args!") return nil
end 


local nthreads = limexDesc.nthreads or 1;

-- create integrator and initialize stages
local limex = LimexTimeIntegrator(nstages)
if ((ndiscs>0 
    and nsolvers>0) 
    or nthreads > 1) then 

  -- multiples discs/solvers 
  if (ndiscs < nstages) then print ("ERROR: Number of discretizations too small:"..ndiscs)  return nil end
  if (nsolvers < nstages) then print ("ERROR: Number of solvers too small:"..nsolvers)  return nil end
  
  for i=1,nstages do 
    limex:add_stage(limexDesc.steps[i], limexDesc.domainDisc[i], limexDesc.nonlinSolver[i])
  end
else 

  -- single disc/solver (=>serial version) 
  for i=1,nstages do 
     limex:add_stage(limexDesc.steps[i], limexDesc.domainDisc, limexDesc.nonlinSolver)
  end
end

-- set tolerance
local tol = limexDesc.tol or 1e-2
limex:set_tolerance(limexDesc.tol)

-- set default time step
limex:set_time_step(limexDesc.dt)    
local dtmin = 1e-4*limexDesc.dt


return limex
end
