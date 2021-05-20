all: multivalued-electrolyte elementary-electrolyte hxy-ecmc hxy-metropolis hxy-gaussian-noise-metropolis xy-ecmc \
		xy-metropolis xy-gaussian-noise-metropolis

elementary-electrolyte:
	cd src/electrolyte/elementary && $(MAKE)

multivalued-electrolyte:
	cd src/electrolyte/multivalued && $(MAKE)

hxy-ecmc:
	cd src/xy_models/hxy/ecmc && $(MAKE)

hxy-metropolis:
	cd src/xy_models/hxy/metropolis/uniform_noise_metropolis && $(MAKE)

hxy-gaussian-noise-metropolis:
	cd src/xy_models/hxy/metropolis/gaussian_noise_metropolis && $(MAKE)

xy-ecmc:
	cd src/xy_models/xy/ecmc && $(MAKE)

xy-metropolis:
	cd src/xy_models/xy/metropolis/uniform_noise_metropolis && $(MAKE)

xy-gaussian-noise-metropolis:
	cd src/xy_models/xy/metropolis/gaussian_noise_metropolis && $(MAKE)
