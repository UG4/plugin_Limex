-- Extend path for packages (if loaded by require)
-- package.path = package.path..";".. ug_get_current_path().."?.lua"
-- print("LIMEX-Plugin: Extended package.path:"..package.path)

-- Extend utilities
util = util or {}
if not(util.limex) then
  ug_load_script("limex_util.lua")
end