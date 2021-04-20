all: multivalued-electrolyte elementary-electrolyte hxy-ecmc hxy-metropolis xy-ecmc xy-metropolis

elementary-electrolyte:
	cd src/electrolyte/elementary && $(MAKE)

multivalued-electrolyte:
	cd src/electrolyte/multivalued && $(MAKE)

hxy-ecmc:
	cd src/xy_models/hxy/ecmc && $(MAKE)

hxy-metropolis:
	cd src/xy_models/hxy/metropolis && $(MAKE)

xy-ecmc:
	cd src/xy_models/xy/ecmc && $(MAKE)

xy-metropolis:
	cd src/xy_models/xy/metropolis && $(MAKE)
