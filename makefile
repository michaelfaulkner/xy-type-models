all: multivalued-electrolyte elementary-electrolyte hxy-ecmc hxy-metropolis xy-ecmc xy-metropolis

multivalued-electrolyte:
	cd src/maggs_electrolyte/multivalued_charges && $(MAKE)

elementary-electrolyte:
	cd src/maggs_electrolyte/elementary_charges && $(MAKE)

hxy-ecmc:
	cd src/hxy/ecmc && $(MAKE)

hxy-metropolis:
	cd src/hxy/metropolis && $(MAKE)

xy-ecmc:
	cd src/xy/ecmc && $(MAKE)

xy-metropolis:
	cd src/xy/metropolis && $(MAKE)
