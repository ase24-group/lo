local about=[[

lib: LUA is a "batteries not included" language. Here are my batteries.
(c) 2023, Tim Menzies, BSD-2. Post issues to http://github.com/timm/lo

USAGE:
  lib=require"lib"]]

local l={}

-- ## Lint

-- Cache things known before.
local b4={}; for k,_ in pairs(_ENV) do b4[k]=k end

-- Complain if we have messed with what was before.
function l.rogues()
  for k,v in pairs(_ENV) do
    if not b4[k] then print("-- W: rogue",k,type(v)) end end end

-- ## Numbers

-- Round to `ndecs` decimals.
function l.rnd(n, ndecs)
  if type(n) ~= "number" then return n end
  if math.floor(n) == n  then return n end
  local mult = 10^(ndecs or 3)
  return math.floor(n * mult + 0.5) / mult end

-- Return a sample from a Gaussian of  `mu` and  `sd`
function l.gaussian(mu,sd)
  mu,sd = mu or 0, sd or 1
  local sq,pi,log,cos,r = math.sqrt,math.pi,math.log,math.cos,math.random
  return mu + sd * sq(-2*log(r())) * cos(2*pi*r())  end

-- Convert list of numbers to mean, standard deviation and sample size. For more on this see
-- [here](http://datagenetics.com/blog/november22017/index.html).
-- This function also returns the mean and sample size of this list 
function l.t2stats(t,    d,n,mu,m2)
	n,mu,m2 = 0,0,0
  for _,x in pairs(t) do
    n  = n+1
    d  = x-mu
    mu = mu + d/n
    m2 = m2 + d*(x-mu) end
	return mu, n<2 and 0 or (m2/(n - 1))^.5, n end

-- ## Lists

-- Return self
function l.self(x) return x end

-- Return any item.
function l.any(t) return t[math.random(#t)] end

-- Return any `n` items (there may be repeats).
function l.many(t,  n,     u)
  n = n or #t
  u={}; for _ = 1,n do u[1+#u] = l.any(t) end; return u end

-- Push onto a list `t`.
function l.push(t,x) t[1+#t]=x; return x end

-- Deep copy
function l.copy(t,    u)
  if type(t) ~= "table" then return t end
  u={}; for k,v in pairs(t) do u[l.copy(k)]=l.copy(v) end
  return u end

-- Lua's default sort does not return the sorted list. So....
function l.sort(t,  fun)
  table.sort(t,fun); return t end

-- Return a function taht sorts on field `x`.
function l.lt(x)
  return function(a,b) return a[x] < b[x] end end

-- Return the `p`-th item in `t`.
function l.per(t, p) return t[(p * #t) // 1] end

-- Return the middle thing.
function l.median(t) return l.per(t, .5) end

-- Return  spread of the data
-- (see 
-- [computing percentiles](https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_probability/bs704_probability10.html)).
function l.stdev(t)   return (l.per(t,.9) - l.per(t,.1))/2.564 end

-- ent is the how mx binary chops are needed to isolate all the parts of a signal.
function l.ent(t, e, n)
    e, n = 0, 0
    for _, v in pairs(t) do n = n + v end
    for _, v in pairs(t) do e = e - v / n * math.log(v / n, 2) end
    return e end

-- Ensure `t` has the fields in `defaults`.
function l.defaults(t, defaults)
  t = t or {}
  for k,v in pairs(defaults) do
    if t[k] == nil then t[k] = v end end
  return t end

-- Apply `fun` to all items in `t`.
function l.map(t,fun,...) --> t
  local u={};  for _,v in pairs(t) do u[1+#u] = fun(v,...) end; return u end

-- Apply `fun` to all keys and values in `t`.
function l.kap(t,fun,...)
  local u = {}; for k, v in pairs(t) do
                  u[1+#u] = fun(k,v,...) end; return u end

-- Return keys and values inside `t` in alphabetical order of the keys.
function l.items(t,fun,    u,i)
  u={}; for k,_ in pairs(t) do u[1+#u]=k end
  table.sort(u,fun)
  i=0
  return function()
    if i<#u then i=i+1; return u[i], t[u[i]] end end end

-- Schwartzian transform:  decorate, sort, undecorate
function l.keysort(t,fun,      tmp)
  tmp = l.map(t,  function(x) return {x=x, y=fun(x)} end) -- decorate
  tmp = l.sort(tmp, l.lt"y")                             -- sort
  return l.map(tmp, function(xy) return xy.x end) end    -- undecorate

-- Given a list of tables all with the same keys, then (a) print the
-- the keys (in alphabetical order); then (b) print the values of
-- each table.
function l.report(ts, header, nwidth, u, say)
    print(header)
    function say(x) io.write(l.fmt("%" .. (nwidth or 4) .. "s", x)) end
    u = {}
    for _, t in pairs(ts) do
      for k, _ in pairs(t) do u[1 + #u] = k end
      table.sort(u)
      say ""; for _, k in pairs(u) do say(k) end; print ""
      for k1, t in l.items(ts) do
        say(k1); for _, k2 in pairs(u) do say(t[k2]) end; print "" end
      return 1 end end

-- ## Thing to string

-- Like printf.
l.fmt = string.format

-- Turn `x` to a string, print it.
function l.oo(x,  ndecs)
  print(l.o(x,ndecs)); return x end

-- Turn a nested list `x`  into a string.
function l.o(x,  ndecs,     fun, u)
  function fun(k, v)
    k = tostring(k)
    if not k:find "^_" then
      return l.fmt(":%s %s", k, l.o(v, ndecs)) end end
  if type(x) == "number" then return tostring(l.rnd(x,ndecs)) end
  if type(x) ~= "table" then return tostring(x) end
  u = #x == 0 and l.sort(l.kap(x, fun)) or l.map(x, l.o, ndecs)
  return "{"..table.concat(u,", ").."}" end

-- ## String to thing

-- String to int or float or nil or bool.
function l.coerce(s1,    fun)
  function fun(s2)
    if s2=="nil" then return nil
    else return s2=="true" or (s2~="false" and s2) end end
  return math.tointeger(s1) or tonumber(s1) or fun(s1:match'^%s*(.*%S)') end

-- String to list of items.
function l.cells(s1,    t)
  t={}; for s2 in s1:gmatch("([^,]+)") do t[1+#t]=l.coerce(s2) end;
  return t end

-- Iterate over a csv file, one row as a time.
function l.csv(src)
  src = src=="-" and io.stdin or io.input(src)
  return function(   line)
    line = io.read()
    if line then return l.cells(line) else io.close(src) end end end

-- For each `key=value` in `t`, look for a command-line flag `-k` or `--key`
-- then update `value` from command line. If the old `value` is a boolean,
-- the `-k` is enough to flip it 
function l.cli(t)
  for k,v in pairs(t) do
    v = tostring(v)
    for n,x in ipairs(arg) do
      if x=="-"..(k:sub(1,1)) or x=="--"..k then
        v= ((v=="false" and "true") or (v=="true" and "false") or arg[n+1])
        t[k] = l.coerce(v) end end end
  if t.help then os.exit(print(t._help)) end
  return t end

-- Turn help string into a options array.
function l.settings(s,    t,pat)
  t={}
  pat = "\n[%s]+[-][%S][%s]+[-][-]([%S]+)[^\n]+= ([%S]+)"
  for k, s1 in s:gmatch(pat) do t[k] = l.coerce(s1) end
  t._help = s
  return t,s end

-- ## Run demos

function l.toplevel ()
  return not pcall(debug.getlocal,5,1) end

-- Run one.
function l.try(s, settings, fun,     before,status)
  math.randomseed(settings.seed or 1234567891)
  before={}; for k,v in pairs(settings) do b4[k]=v end
  io.write("🔷 ".. s.." ")
  status = fun()==false
  for k,v in pairs(before) do settings[k]=v end
  if   status
  then print(" ❌ FAIL"); return true
  else print("✅ PASS"); return false end  end

-- Run all the requests on the command line.
function l.run(settings, funs)
  l.cli(settings)
  for _,com in pairs(arg) do
     if com=="all" then  l.runall(settings,funs) end
     if funs[com] then l.try(com, settings, funs[com]) end end
  l.rogues() end

-- Run all.
function l.runall(settings,funs,     oops)
  oops = -1 -- we have one test that deliberately fails
  for k,fun in l.items(funs) do
    if k~="all" then
      if l.try(k,settings, fun) then oops = oops + 1 end end end
  l.rogues()
  os.exit(oops) end

-- ##  Discrete Performance Stats
-- NEW

-- ### ABCD
-- For one class, calcuate statistics for symbolic classification.

-- Create,
function l.ABCD(klass, before)
  return {klass=klass, a=(before or 0), b=0, c=0, d=0} end

-- Update.
function l.abcd(abcd1, want,got)
  if   want == abcd1.klass
  then if want==got       then abcd1.d=abcd1.d+1 else abcd1.b=abcd1.b+1 end
  else if got==abcd1.klass then abcd1.c=abcd1.c+1 else abcd1.a=abcd1.a+1 end end end

-- Query.
function l.pf(abcd1)        return abcd1.c           / (abcd1.a+abcd1.c+1E-30) end
function l.recall(abcd1)    return abcd1.d           / (abcd1.b+abcd1.d+1E-30) end
function l.accuracy(abcd1)  return (abcd1.a+abcd1.d) / (abcd1.a+abcd1.b+abcd1.c+abcd1.d+1E-30) end
function l.precision(abcd1) return abcd1.d           / (abcd1.c+abcd1.d+1E-30) end
function l.f(abcd1,   p,r)  p,r  = l.precision(abcd1),l.recall(abcd1); return (2*p*r)  / (p+r) end
function l.g(abcd1,   nf,r)  nf,r = 1-l.pf(abcd1),l.recall(abcd1);     return (2*nf*r) / (nf+r) end

-- ### ABCDS
-- For many classes,  calcuate statistics for symbolic classification.

-- Create.
function l.ABCDS() return {all={}, n=0} end

-- Update.
function l.abcds(abcds1,want,got)
  abcds1.all[want] = abcds1.all[want] or l.ABCD(want, abcds1.n)
  abcds1.n         = abcds1.n + 1
  for _,abcd1 in pairs(abcds1.all) do l.abcd(abcd1,want,got) end end

-- Report.
function l.abcdsreport(abcds1,     u)
  u={}; for k,abcd1 in pairs(abcds1.all) do u[k] = {
    _n   = abcd1.a+abcd1.b+abcd1.c+abcd1.d,
    _a   = abcd1.a,
    _b   = abcd1.b,
    _c   = abcd1.c,
    _d   = abcd1.d,
    acc  = l.accuracy(abcd1),
    prec = l.precision(abcd1),
    pd   = l.recall(abcd1),
    pf   = l.pf(abcd1),
    f    = l.f(abcd1),
    g    = l.g(abcd1)} end
  return u end

-- Non-parametric effect-size test
-- M.Hess, J.Kromrey. 
-- Robust Confidence Intervals for Effect Sizes: 
-- A Comparative Study of Cohen's d and Cliff's Delta Under Non-normality and Heterogeneous Variances
-- American Educational Research Association, San Diego, April 12 - 16, 2004    
-- 0.147=  small, 0.33 =  medium, 0.474 = large; med --> small at .2385
function l.cliffsDelta(ns1,ns2,  cliffs) 
  if #ns1 > 256     then ns1 = l.many(ns1,256) end
  if #ns2 > 256     then ns2 = l.many(ns2,256) end
  if #ns1 > 10*#ns2 then ns1 = l.many(ns1,10*#ns2) end
  if #ns2 > 10*#ns1 then ns2 = l.many(ns2,10*#ns1) end
  local n,gt,lt = 0,0,0
  for _,x in pairs(ns1) do
    for _,y in pairs(ns2) do
      n = n + 1
      if x > y then gt = gt + 1 end
      if x < y then lt = lt + 1 end end end
  return math.abs(lt - gt)/n > cliffs end

-- Non-parametric significant test. From p220 to 223 of theg
-- Efron text  'introduction to the boostrap'.
-- https://bit.ly/3iSJz8B Typically, conf=0.05 and bootstraps is 100s to 1000s.
function l.bootstrap(y0,z0,bootstraps,conf,      _,__)
  local function delta(y,z)
    local ymu,ysd,yn = l.t2stats(y)
    local zmu,zsd,zn = l.t2stats(z)
    return math.abs(ymu - zmu) / ((1E-30 + ysd^2/yn + zsd^2/zn)^.5)
  end ------------------------------------------------------------ 
  local x, y, z, yhat, zhat = l.copy(y0), l.copy(y0), l.copy(z0), {}, {}
  for _,z1 in pairs(z0) do l.push(x,z1) end
  local n, tobs = 0, delta(y,z)
  local xmu,ymu,zmu = l.t2stats(x), l.t2stats(y), l.t2stats(z)
  for _,y1 in pairs(y0) do l.push(yhat, y1 - ymu + xmu) end
  for _,z1 in pairs(z0) do l.push(zhat, z1 - zmu + xmu) end
  for _= 1,bootstraps do
    if delta(l.many(yhat), l.many(zhat)) > tobs then n = n + 1 end end 
  return n / bootstraps >= conf end

  -- ## Main
return l
