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

local limexDesc = {

  nstages = 2,
  steps = {1,2,3},
  domainDisc={domainDisc[1], domainDisc[2], domainDisc[3]}
  nonlinSolver = {nlSolver[1], nlSolver[2], nlSolver[3]}
 --lSolver = limexLSolver,
  
  tol = tol,
  dt = dtlimex,
  dtmin = 1e-9,
  
}

--]]

ug_load_script("util/table_desc_util.lua")
ug_load_script("util/table_util.lua")

util = util or {}
util.limex = util.limex or {}


-- function for creating an integrator
function util.limex.CreateIntegrator(limexDesc)

local nstages = limexDesc.nstages
local nsteps = table.getn(limexDesc.steps)
local ndiscs = table.getn(limexDesc.domainDisc)
local nsolvers = table.getn(limexDesc.nonlinSolver)

if ((nsteps < nstages) or (ndiscs < nstages) or (nsolvers < nstages)) then 
  print ("ERROR: Array too short!")
  return nil
end

-- create integrator and push stages
local limex = LimexTimeIntegrator(nstages)
for i=1,nstages do 
  limex:add_stage(limexDesc.steps[i], limexDesc.domainDisc[i], limexDesc.nonlinSolver[i])
end

-- tolerance
local tol = limexDesc.tol or 1e-2
limex:set_tolerance(limexDesc.tol)

-- default time step
limex:set_time_step(limexDesc.dt)    
local dtmin = 1e-4*limexDesc.dt


return limex
end
