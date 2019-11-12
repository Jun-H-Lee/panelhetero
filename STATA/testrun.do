/// test run for Stata package panelhetero

clear 

use panelinit

xtset id time

/// T = 50

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")

/// T = 49
drop if time == 50

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")

/// T = 48
drop if time == 49

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")

/// T = 47
drop if time == 48

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")

/// T = 46
drop if time == 47

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")

/// T = 45
drop if time == 46

phkd data, method("naive")
phkd data, method("hpj")
phkd data, method("toj")

phmoment data, method("naive")
phmoment data, method("hpj")
phmoment data, method("toj")

phecdf data, method("naive")
phecdf data, method("hpj")
phecdf data, method("toj")
