run-all-make-files:
	cd src/2d_maggs_electrolyte/multivalued_charges && $(MAKE)
	cd src/2d_maggs_electrolyte/elementary_charges && $(MAKE)
	cd src/hxy/ecmc && $(MAKE)
	cd src/hxy/metropolis && $(MAKE)
	cd src/xy/ecmc && $(MAKE)
	cd src/xy/metropolis && $(MAKE)
