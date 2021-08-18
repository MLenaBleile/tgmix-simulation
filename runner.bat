	@echo off
	SET SIMDIR=C:\Users\s198663\Documents\TGM\simulation

:repeat
C:/PROGRA~1/R/R-40~1.1/bin/R.exe CMD BATCH tgmix-simulation.R || goto :repeat
echo Success!
   